<!--
CATEGORY: architecture (explanation of internals / how a system is built)

Exemplars in this file and their provenance:
  - Architecture of SQLite            - sqlite.org                     - Public domain (no obligation)
  - rustc: Overview of the compiler   - rust-lang/rustc-dev-guide      - MIT OR Apache-2.0 (attribute)
  - rustc: The compiler source code   - rust-lang/rustc-dev-guide      - MIT OR Apache-2.0 (attribute)
-->

# Architecture

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. Lead with one orienting mental model before any detail (SQLite: "compiles SQL
   into bytecode, then runs that bytecode on a virtual machine"; rustc: the
   compiler is "not a series of sequential passes" but "queries that call each
   other"). One sentence the reader keeps.
2. Decompose into named components and say what each is FOR, then where it lives
   (SQLite maps each subsystem to its source file; rustc maps each stage to its
   rustc_* crate). Structure follows function.
3. Explain design choices by the goal they serve, not just what they are (rustc:
   the query system exists TO enable incremental compilation - redo only what
   changed; SQLite: why a bytecode VM rather than a tree-walker).
4. Surface the hard constraint or tradeoff that shaped the design (rustc: parallelism
   is "generally a hard problem", so most of the compiler is still serial and the
   team is converting RefCells to Mutexes incrementally).
Anti-patterns: a file-by-file inventory with no mental model; describing WHAT
each part is without WHY it exists or how the pieces interact.
-->

## Architecture of SQLite

<!-- Annotation: Famously compact. One-sentence overview (SQL -> bytecode -> VM), then
each component gets a role and a source file. Note the rationale asides - why
Lemon over YACC, why the tokenizer drives the parser - which turn an inventory
into an explanation. -->

```text
Introduction

This document describes the architecture of the SQLite library. The information here is useful to those who want to understand or modify the inner workings of SQLite. 

A nearby diagram shows the main components of SQLite and how they interoperate. The text below explains the roles of the various components. 

Overview

SQLite works by compiling SQL text into bytecode, then running that bytecode using a virtual machine.

The sqlite3_prepare_v2() and related interfaces act as a compiler for converting SQL text into bytecode. The sqlite3_stmt object is a container for a single bytecode program that implements a single SQL statement. The sqlite3_step() interface passes a bytecode program into the virtual machine, and runs the program until it either completes, or forms a row of result to be returned, or hits a fatal error, or is interrupted. 

Interface

Much of the C-language Interface is found in source files main.c, legacy.c, and vdbeapi.c though some routines are scattered about in other files where they can have access to data structures with file scope. The sqlite3_get_table() routine is implemented in table.c. The sqlite3_mprintf() routine is found in printf.c. The sqlite3_complete() interface is in complete.c. The TCL Interface is implemented by tclsqlite.c.

To avoid name collisions, all external symbols in the SQLite library begin with the prefix sqlite3. Those symbols that are intended for external use (in other words, those symbols which form the API for SQLite) add an underscore, and thus begin with sqlite3_. Extension APIs sometimes add the extension name prior to the underscore; for example: sqlite3rbu_ or sqlite3session_.

Tokenizer

When a string containing SQL statements is to be evaluated it is first sent to the tokenizer. The tokenizer breaks the SQL text into tokens and hands those tokens one by one to the parser. The tokenizer is hand-coded in the file tokenize.c.

Note that in this design, the tokenizer calls the parser. People who are familiar with YACC and BISON may be accustomed to doing things the other way around - having the parser call the tokenizer. Having the tokenizer call the parser is better, though, because it can be made threadsafe and it runs faster.

Parser

The parser assigns meaning to tokens based on their context. The parser for SQLite is generated using the Lemon parser generator. Lemon does the same job as YACC/BISON, but it uses a different input syntax which is less error-prone. Lemon also generates a parser which is reentrant and thread-safe. And Lemon defines the concept of a non-terminal destructor so that it does not leak memory when syntax errors are encountered. The grammar file that drives Lemon and that defines the SQL language that SQLite understands is found in parse.y.

Because Lemon is a program not normally found on development machines, the complete source code to Lemon (just one C file) is included in the SQLite distribution in the "tool" subdirectory.
Code Generator

After the parser assembles tokens into a parse tree, the code generator runs to analyze the parse tree and generate bytecode that performs the work of the SQL statement. The prepared statement object is a container for this bytecode. There are many files in the code generator, including: attach.c, auth.c, build.c, delete.c, expr.c, insert.c, pragma.c, select.c, trigger.c, update.c, vacuum.c, where.c, wherecode.c, and whereexpr.c. In these files is where most of the serious magic happens. expr.c handles code generation for expressions. where*.c handles code generation for WHERE clauses on SELECT, UPDATE and DELETE statements. The files attach.c, delete.c, insert.c, select.c, trigger.c update.c, and vacuum.c handle the code generation for SQL statements with the same names. (Each of these files calls routines in expr.c and where.c as necessary.) All other SQL statements are coded out of build.c. The auth.c file implements the functionality of sqlite3_set_authorizer().

The code generator, and especially the logic in where*.c and in select.c, is sometimes called the query planner. For any particular SQL statement, there might be hundreds, thousands, or millions of different algorithms to compute the answer. The query planner is an AI that strives to select the best algorithm from these millions of choices.
Bytecode Engine

The bytecode program created by the code generator is run by a virtual machine.

The virtual machine itself is entirely contained in a single source file vdbe.c. The vdbe.h header file defines an interface between the virtual machine and the rest of the SQLite library and vdbeInt.h which defines structures and interfaces that are private to the virtual machine itself. Various other vdbe*.c files are helpers to the virtual machine. The vdbeaux.c file contains utilities used by the virtual machine and interface modules used by the rest of the library to construct VM programs. The vdbeapi.c file contains external interfaces to the virtual machine such as the sqlite3_bind_int() and sqlite3_step(). Individual values (strings, integer, floating point numbers, and BLOBs) are stored in an internal object named "Mem" which is implemented by vdbemem.c.

SQLite implements SQL functions using callbacks to C-language routines. Even the built-in SQL functions are implemented this way. Most of the built-in SQL functions (ex: abs(), count(), substr(), and so forth) can be found in the func.c source file. Date and time conversion functions are found in date.c. Some functions such as coalesce() and typeof() are implemented as bytecode directly by the code generator.
B-Tree

An SQLite database is maintained on disk using a B-tree implementation found in the btree.c source file. Separate B-trees are used for each table and each index in the database. All B-trees are stored in the same disk file. The file format details are stable and well-defined and are guaranteed to be compatible moving forward.

The interface to the B-tree subsystem and the rest of the SQLite library is defined by the header file btree.h.
Page Cache

The B-tree module requests information from the disk in fixed-size pages. The default page_size is 4096 bytes but can be any power of two between 512 and 65536 bytes. The page cache is responsible for reading, writing, and caching these pages. The page cache also provides the rollback and atomic commit abstraction and takes care of locking of the database file. The B-tree driver requests particular pages from the page cache and notifies the page cache when it wants to modify pages or commit or rollback changes. The page cache handles all the messy details of making sure the requests are handled quickly, safely, and efficiently.

The primary page cache implementation is in the pager.c file. WAL mode logic is in the separate wal.c file. In-memory caching is implemented by the pcache.c and pcache1.c files. The interface between the page cache subsystem and the rest of SQLite is defined by the header file pager.h.
OS Interface

In order to provide portability across operating systems, SQLite uses an abstract object called the VFS. Each VFS provides methods for opening, reading, writing, and closing files on disk, and for other OS-specific tasks such as finding the current time, or obtaining randomness to initialize the built-in pseudo-random number generator. SQLite currently provides VFSes for unix (in the os_unix.c file) and Windows (in the os_win.c file).
Utilities

Memory allocation, caseless string comparison routines, portable text-to-number conversion routines, and other utilities are located in util.c. Symbol tables used by the parser are maintained by hash tables found in hash.c. The utf.c source file contains Unicode conversion subroutines. SQLite has its own private implementation of printf() (with some extensions) in printf.c and its own pseudo-random number generator (PRNG) in random.c.
Test Code

Files in the "src/" folder of the source tree whose names begin with test are for testing only and are not included in a standard build of the library.
```

## rustc - Overview of the compiler (the "big design choice" lens)

<!-- Annotation: The exemplary move is leading with the ONE design decision that
shapes everything else - rustc is organized as demand-driven queries, not
sequential passes - and immediately giving the reason (incremental compilation:
on recompile, redo as little as possible). Only after that mental model is in
place does it walk the pipeline (lexer -> parser -> ... -> codegen). Note also
the honesty about limits: it states plainly that most of the compiler is not yet
parallel because parallelism "is generally a hard problem," which tells the
reader the shape of the codebase as it really is, not an idealized version. -->

```text
The first big implementation choice is Rust's use of the query system in its
compiler. The Rust compiler is not organized as a series of passes over the
code which execute sequentially. The Rust compiler does this to make
incremental compilation possible -- that is, if the user makes a change to
their program and recompiles, we want to do as little redundant work as
possible to output the new binary.

In rustc, all the major steps above are organized as a bunch of queries that
call each other. For example, there is a query to ask for the type of
something and another to ask for the optimized MIR of a function. These
queries can call each other and are all tracked through the query system, and
the results of the queries are cached on disk so that the compiler can tell
which queries' results changed from the last compilation and only redo those.
This is how incremental compilation works.

[Pipeline, in brief]
Command line argument parsing occurs in the rustc_driver. This crate defines
the compile configuration that is requested by the user and passes it to the
rest of the compilation process as a rustc_interface::Config.

The raw Rust source text is analyzed by a low-level lexer located in
rustc_lexer. At this stage, the source text is turned into a stream of atomic
source code units known as tokens. The token stream passes through a
higher-level lexer located in rustc_parse to prepare for the next stage of the
compile process. ... The rustc_ast::ast::{Crate, Expr, Pat, ...} AST nodes are
returned from the parser while the standard Diag API is used for error
handling. Generally Rust's compiler will try to recover from errors by parsing
a superset of Rust's grammar, while also emitting an error type.

We then begin what is simply called code generation or codegen. ... Since rustc
uses LLVM for code generation, the first step is to convert the MIR to LLVM-IR.
This is where the MIR is actually monomorphized. The LLVM-IR is passed to LLVM,
which does a lot more optimizations on it, emitting machine code.

[On limits, stated honestly]
Currently, there is only one part of rustc that is parallel by default: code
generation. However, the rest of the compiler is still not yet parallel. There
have been lots of efforts spent on this, but it is generally a hard problem.
```

## rustc - The compiler source code (the "repository structure" lens)

<!-- Annotation: A different architectural view of the same system - not the data
flow but the dependency STRUCTURE. Note how it explains the shape as a
consequence of a goal: the query system is wired with function pointers
specifically to break inter-crate dependencies and enable more parallel
compilation. It also surfaces the real-world cost of a design (rustc_middle is
"a really large crate, leading to long compile times") rather than presenting the
structure as frictionless. Architecture writing that admits tradeoffs is more
useful than architecture writing that only celebrates. -->

```text
At the very bottom of the dependency tree are a handful of crates that are used
by the whole compiler (e.g. rustc_span). The very early parts of the
compilation process (e.g. parsing and the AST) depend on only these.

After the AST is constructed and other early analysis is done, the compiler's
query system gets set up. The query system is set up in a clever way using
function pointers. This allows us to break dependencies between crates,
allowing more parallel compilation. The query system is defined in rustc_middle,
so nearly all subsequent parts of the compiler depend on this crate. It is a
really large crate, leading to long compile times. Some efforts have been made
to move stuff out of it with varying success.

At the top of the dependency tree is rustc_driver and rustc_interface which is
an unstable wrapper around the query system helping drive various stages of
compilation. Other consumers of the compiler may use this interface in
different ways (e.g. rustdoc or maybe eventually rust-analyzer).

The rust-lang/rust repository consists of a single large cargo workspace
containing the compiler, the standard libraries, and rustdoc, along with the
build system and a bunch of tools. The repository consists of three main
directories: compiler/ contains the source code for rustc (around 50
interdependent crates, all named rustc_*); library/ contains the standard
libraries (core, alloc, std, proc_macro, test); and the build system and tools
live alongside them.
```
