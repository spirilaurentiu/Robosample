<!--
CATEGORY: documentation (tutorials, how-tos, and conceptual explanation - teaching
a reader to DO or UNDERSTAND something, as distinct from reference lookup)

Exemplars in this file and their provenance:
  - SQLite "Datatypes / Flexible typing is a feature" - sqlite.org                - Public domain
  - Rust By Example: Vectors                          - rust-lang/rust-by-example - MIT OR Apache-2.0 (attribute)
Attribution for the Rust material: "Rust By Example", rust-lang, dual MIT/Apache-2.0.

WHAT THIS CATEGORY IS FOR vs reference.md - read before using:
  documentation = LEARNING / TASK / UNDERSTANDING. Read start-to-finish. Answers
    "how do I accomplish X?" (how-to), "teach me to use Y" (tutorial), and
    "why does this work this way?" (explanation).
  reference     = LOOKUP. Read by jumping to one entry. Answers "what exactly does
    this function take, return, and error with?"
  Both may mention the SAME function. The difference is the reader's goal:
    - "what does Vec::push return / when does it panic?"      -> reference voice
    - "how do I build up a vector and iterate it?"            -> documentation voice (tutorial)
    - "why does SQLite accept a string in an INT column?"     -> documentation voice (explanation)
-->

# Documentation (tutorials, how-tos & explanation)

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. Orient to the reader's goal first: what they'll be able to do/understand, and
   any prerequisites. A tutorial has a destination.
2. One concept at a time, in dependency order. Never use a thing before showing it.
3. Anticipate the misconception at the exact point it arises and defuse it inline
   (SQLite names the reader's "this must be a bug!" reaction and reframes it).
4. Every actionable claim gets a runnable example, and the example's OUTPUT is
   shown so the reader can check themselves (Rust By Example prints results and
   even comments the lines that intentionally fail).
5. Motivate before mechanism: say why a step matters before the how.
6. Progressive disclosure: simple path first; edge cases and options after,
   clearly marked as such.
Anti-patterns: dumping the full API surface flat (that's reference, not a
tutorial); examples with no expected output; introducing a term before defining
it; a wall of prose with nowhere for the reader to act.
-->

## SQLite - "Flexible typing is a feature, not a bug" (explanation, public domain)

<!-- Annotation: The exemplary move is meeting the reader's WRONG mental model
head-on. A developer coming from Postgres sees a string accepted into an INT
column and concludes SQLite is broken. Instead of burying that, the doc states
the reader's reaction ("a bug") and reframes it as a deliberate design ("a
feature"), explains the actual rule (type affinity is *recommended, not
required*), grounds it in a concrete CREATE/INSERT example, and even points to
the escape hatch for readers who disagree (STRICT tables). Explanation-genre
writing changes how the reader THINKS about a thing; it earns that by naming the
misconception first. -->

```text
Most SQL database engines use static, rigid typing: the datatype of a value is
determined by its container - the column it is stored in. SQLite uses a more
general dynamic type system: the datatype of a value is associated with the
value itself, not with its container. Flexible typing is a feature of SQLite,
not a bug.

SQL engines with rigid typing try to convert values to the column's type.
Consider:

    CREATE TABLE t1(a INT, b VARCHAR(10));
    INSERT INTO t1(a,b) VALUES('123', 456);

A rigidly-typed database converts the string '123' to the integer 123 and the
integer 456 to the string '456' before inserting. So that this example works in
SQLite the same way it works elsewhere, SQLite supports "type affinity" on
columns. The type affinity of a column is the *recommended* type for data
stored there. The important idea: the type is recommended, not required. Any
column can still store any type of data; some columns, given the choice, simply
prefer one storage class over another.

(As an accommodation for those who expect rigid typing, SQLite 3.37.0
introduced STRICT tables, which impose mandatory datatype constraints.)
```

## Rust By Example - Vectors (tutorial, MIT/Apache)

<!-- Annotation: A tutorial doing the "show, run, verify" loop. It introduces one
operation at a time (collect -> vec! -> push -> len -> index -> iterate), and each
line is runnable with its result visible. Note the teaching use of DELIBERATE
failures: the lines that would panic (pushing to an immutable vec; indexing out
of bounds) are shown and commented "FIXME ^ Comment out this line," so the reader
learns the boundary by seeing it, not by being told about it in the abstract.
Contrast this with reference.md's Vec::push entry, which states the panic
condition as a contract - same function, different genre, different job. -->

```text
// Iterators can be collected into vectors
let collected_iterator: Vec<i32> = (0..10).collect();
println!("Collected (0..10) into: {:?}", collected_iterator);
//=> Collected (0..10) into: [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]

// The `vec!` macro can be used to initialize a vector
let mut xs = vec![1i32, 2, 3];
println!("Initial vector: {:?}", xs);
//=> Initial vector: [1, 2, 3]

// Insert a new element at the end of the vector
xs.push(4);
println!("Vector: {:?}", xs);
//=> Vector: [1, 2, 3, 4]

// The `len` method yields the number of elements currently stored
println!("Vector length: {}", xs.len());     //=> 4

// Indexing is done with square brackets (indexing starts at 0)
println!("Second element: {}", xs[1]);       //=> 2

// Out-of-bounds indexing yields a panic - shown on purpose so you see the edge:
// println!("Fourth element: {}", xs[3]);
//   ^ Comment this in and it panics: index out of bounds.

// Vectors can be iterated over
for x in xs.iter() {
    println!("> {}", x);
}
```

---

## Negative pair (bad -> good)

<!-- Original, written for this corpus. The AI-tutorial failure mode: a flat
feature/API dump with no destination, no expected output, and no anticipation of
the one thing that always trips people up. Note the bad version is really
reference material wearing a tutorial's hat - which is exactly the doc-vs-reference
confusion this file's header warns about. -->

**BAD** - flat API dump, no goal, no output, misconception ignored:

```text
## Using the Config Loader

The config loader loads config. It supports JSON, YAML, and env vars. You can
call load(), merge(), and validate(). It also supports overrides and defaults
and nested keys and arrays. Here is the API:

  load(path)
  merge(a, b)
  validate(schema)
  get(key)
  set(key, value)
```

**GOOD** - goal first, one step at a time, output shown, gotcha anticipated:

```text
## Loading your first config

By the end of this you'll load a YAML config, override one value from an
environment variable, and confirm the result. You'll need the package installed
and a file to load.

Start with a minimal config.yaml:

    port: 8080
    logLevel: info

Load it:

    const cfg = await load("config.yaml");
    console.log(cfg.get("port"));   //=> 8080

Now override the log level without editing the file - set an env var and reload:

    LOG_LEVEL=debug node app.js
    console.log(cfg.get("logLevel"));   //=> debug

A common surprise here: env overrides win over file values, not the other way
around. People expect the file they just edited to take precedence. It doesn't -
the precedence order is defaults < file < environment, so deploys can override
committed config without a code change. Keep that order in mind when a value
isn't what you expect.

Next: layering multiple files with merge() ->
```

<!-- Why the good one works: it opens with what the reader will be able to do and
what they need; introduces exactly one capability at a time and shows the output
after each step so the reader can verify; and stops at the precise point of
predictable confusion (env-vs-file precedence) to name the wrong expectation and
explain the real rule with its rationale - the same "name the misconception"
move SQLite makes. The bad version lists the whole API surface, which belongs in
reference.md, not a tutorial. -->