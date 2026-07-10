<!--
CATEGORY: readme (project overview / first-contact document)

Exemplars in this file and their provenance:
  - SQLite "About" + "Distinctive Features" - sqlite.org         - Public domain
  - ripgrep README (excerpt)                - BurntSushi/ripgrep - MIT OR Unlicense
-->

# READMEs / project overviews

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. Answer "what is this?" in the first sentence, in plain terms, before any
   install command or badge (SQLite: "an in-process library that implements a
   self-contained, serverless, zero-configuration SQL database engine").
2. Say who it's for and why they'd choose it - the differentiators, honestly
   (ripgrep leads with recursive-search-by-default and gitignore-aware filtering).
3. Include a "why NOT this" / limitations section. Mature READMEs tell you when
   to reach for something else (ripgrep: "you need a portable, POSIX tool -> grep").
4. Show, then tell: one concrete usage example early beats a paragraph of prose.
5. Order by reader need: what-it-is -> why-care -> quickstart -> depth. Not
   contributor/build details first.
Anti-patterns: opening with badges and TOC before saying what the thing is;
feature lists with no differentiation; no limitations section; burying the
one-line description under installation instructions.
-->

## SQLite - "About" (the what-it-is, done well)

<!-- Annotation: The first sentence is a complete, jargon-controlled definition.
Every clause earns its place: "in-process", "self-contained", "serverless",
"zero-configuration", "transactional" each rule out a category of alternative.
Then it immediately grounds the abstraction ("reads and writes directly to
ordinary disk files") and states the payoff (single-file, cross-platform). No
badges, no preamble - the reader knows what SQLite is by the end of line one. -->

```text
SQLite is an in-process library that implements a self-contained,
serverless, zero-configuration, transactional SQL database engine. The
code for SQLite is in the public domain and is thus free for use for any
purpose, commercial or private.

Unlike most other SQL databases, SQLite does not have a separate server
process. SQLite reads and writes directly to ordinary disk files. A
complete SQL database with multiple tables, indices, triggers, and views,
is contained in a single disk file. The database file format is
cross-platform - you can freely copy a database between 32-bit and 64-bit
systems or between big-endian and little-endian architectures. These
features make SQLite a popular choice as an Application File Format.

SQLite is a compact library. With all features enabled, the library size
can be less than 900KiB. There is a tradeoff between memory usage and
speed. SQLite generally runs faster the more memory you give it.
Nevertheless, performance is usually quite good even in low-memory
environments.
```

## SQLite - "Distinctive Features" (the why-choose-this)

<!-- Annotation: Differentiation by contrast. Instead of "SQLite is easy", it
shows the thing competitors make you do and that SQLite doesn't: no install, no
server process, no config files, nothing to recover after a crash. The genre
move: define yourself against the reader's existing mental model. -->

```text
SQLite does not need to be "installed" before it is used. There is no
"setup" procedure. There is no server process that needs to be started,
stopped, or configured. There is no need for an administrator to create a
new database instance or assign access permissions to users. SQLite uses
no configuration files. Nothing needs to be done to tell the system that
SQLite is running. No actions are required to recover after a system crash
or power failure. There is nothing to troubleshoot. SQLite just works.

Other more familiar database engines run great once you get them going.
But doing the initial installation and configuration can be intimidatingly
complex.
```

## ripgrep - README (differentiators + honest "why not")

<!-- Annotation: Two things to copy here. First, the differentiators are specific
and comparative, not adjectival ("defaults to recursive search", "won't search
files ignored by your .gitignore", "won't search hidden or binary files") - each
is a concrete behavior. Second, the rare and valuable "Why shouldn't I use
ripgrep?" section: a mature README that names its own limits earns trust and
saves the reader time. -->

```text
ripgrep is a line-oriented search tool that recursively searches the
current directory for a regex pattern. By default, ripgrep will respect
gitignore rules and automatically skip hidden files/directories and binary
files. ripgrep has first class support on Windows, macOS and Linux, with
binary downloads available for every release.

Like other tools specialized to code search, ripgrep defaults to recursive
search and does automatic filtering. Namely, ripgrep won't search files
ignored by your .gitignore/.ignore/.rgignore files, it won't search hidden
files and it won't search binary files. Automatic filtering can be disabled
with rg -uuu.

Why shouldn't I use ripgrep? ... The primary reasons not to use ripgrep
probably consist of one or more of the following:

  - You need a portable and ubiquitous tool. While ripgrep works on
    Windows, macOS and Linux, it is not ubiquitous and it does not conform
    to any standard such as POSIX. The best tool for this job is good old
    grep.
  - There still exists some other feature (or bug) not listed in this
    README that you rely on that's in another tool that isn't in ripgrep.
  - There is a performance edge case where ripgrep doesn't do well where
    another tool does do well.
```

## Negative pair (bad -> good)

<!-- Original, written for this corpus - no external source. Teaches the boundary
directly: the failure mode of AI-generated READMEs is starting with scaffolding
and adjectives instead of a concrete definition and differentiators. -->

**BAD** - badges-and-adjectives opening, no real definition:

```text
# FooCache 🚀

![build](badge) ![coverage](badge) ![license](badge)

## Table of Contents
1. Installation
2. Usage
3. Contributing

FooCache is a powerful, modern, blazing-fast, easy-to-use caching solution
built with cutting-edge technology to supercharge your application's
performance. Get started in seconds!

## Installation
npm install foocache
```

**GOOD** - definition first, differentiator, honest scope:

```text
# FooCache

FooCache is an in-process LRU cache for Node.js with a bounded memory
footprint: you set a byte budget, and it evicts by true serialized size
rather than entry count. It's for services that cache heterogeneous objects
(a few bytes to a few MB) where counting entries doesn't bound memory.

If your entries are uniform in size, a simpler count-based cache (e.g.
`lru-cache`) is lighter and probably a better fit.

    import { FooCache } from "foocache";
    const cache = new FooCache({ maxBytes: 64 * 1024 * 1024 });
    cache.set("user:42", user);   // evicts oldest when the 64 MB budget is hit

## Installation
    npm install foocache
```

<!-- Why the good one works: sentence one is a real definition with a bound that
rules out alternatives; it names the specific reader it's for; it points elsewhere
when FooCache is the wrong choice; and it shows the distinguishing behavior in
five lines before any install boilerplate. -->