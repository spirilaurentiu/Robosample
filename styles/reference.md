<!--
CATEGORY: reference (lookup material - you arrive knowing what you want, you leave
with the exact signature/behavior/return codes/errors/edge cases)

Exemplars in this file and their provenance:
  - SQLite C API: sqlite3_step()      - sqlite.org                 - Public domain (no obligation)
  - Rust std: Vec / Vec::push         - rust-lang/rust (std docs)  - MIT OR Apache-2.0 (attribute)
Attribution for the Rust material: "The Rust Standard Library documentation",
rust-lang, dual MIT/Apache-2.0.
-->

# Reference documentation

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. Lead with the exact signature/prototype. No throat-clearing - show
   `int sqlite3_step(sqlite3_stmt*)` first, then explain.
2. Enumerate EVERY return/error case and what each means, including the
   annoying ones (SQLite documents SQLITE_BUSY/DONE/ROW/ERROR/MISUSE each with
   its own recovery action).
3. State edge conditions, preconditions, and complexity precisely: what must be
   called first, what happens on misuse, and cost (Rust: push is "amortized O(1)";
   "Panics if the new capacity exceeds isize::MAX bytes").
4. Be honest about warts. SQLite literally labels a legacy misfeature "Goofy
   Interface Alert" and explains the trap - reference docs earn trust by naming
   sharp edges instead of hiding them.
5. One tight example per entry showing the common call and its result, not a
   tutorial (Rust: `let mut vec = vec![1, 2]; vec.push(3); assert_eq!(vec, [1,2,3]);`).
6. Uniform structure across entries so the reader learns the shape once and scans.
Anti-patterns: motivational prose instead of the contract; "usually returns X"
where the contract should be exact; missing error/panic cases; examples that
teach a workflow instead of demonstrating the one call.
-->

## SQLite C API - `sqlite3_step()` (public domain)

<!-- Annotation: A model reference entry. Signature first; then the return set is
enumerated exhaustively, and crucially each code comes with the caller's required
response (BUSY -> retry or rollback; DONE -> reset before reusing; ROW -> read the
columns; MISUSE -> your code is broken). Note the "Goofy Interface Alert": rather
than hide a legacy design wart, it flags it, explains the trap (generic
SQLITE_ERROR until you call reset/finalize), and points at the fix (the v2
interface). Documenting the sharp edge IS the service. -->

```text
int sqlite3_step(sqlite3_stmt*);

After a prepared statement has been prepared using sqlite3_prepare_v2() (or a
sibling), this function must be called one or more times to evaluate the
statement.

Return codes (v2 interface - the "vX" prepares return these directly):

  SQLITE_ROW    A new row of result data is ready. Read it with the column
                access functions, then call sqlite3_step() again for the next
                row.
  SQLITE_DONE   The statement has finished executing successfully.
                sqlite3_step() should not be called again without first calling
                sqlite3_reset() to reset the virtual machine to its initial
                state.
  SQLITE_BUSY   The engine could not acquire the locks it needs. If the
                statement is a COMMIT or occurs outside an explicit
                transaction, you can retry it. If it is not a COMMIT and occurs
                inside an explicit transaction, you should ROLLBACK and retry.
  SQLITE_ERROR  A run-time error (e.g. a constraint violation). Do not call
                sqlite3_step() again on this statement; call sqlite3_errmsg()
                for details.
  SQLITE_MISUSE The routine was called inappropriately - e.g. on a statement
                that was already finalized, or from two threads on one
                connection at once. An application that receives SQLITE_MISUSE
                is broken by definition.

Goofy Interface Alert: In the legacy interface, sqlite3_step() always returns a
generic SQLITE_ERROR following any error other than SQLITE_BUSY and
SQLITE_MISUSE; you must call sqlite3_reset() or sqlite3_finalize() to recover
the specific error code. This is a goofy design. It is fixed in the v2
interface, which returns the specific code directly - use the v2 prepares.
```

## Rust std - `Vec` / `Vec::push` (MIT/Apache)

<!-- Annotation: Note the section discipline the Rust std docs apply to every
method: a one-line summary of what it does, an explicit "Panics" clause stating
exactly when it aborts, a complexity note, and a runnable "Examples" block whose
asserts double as a spec. The surrounding type docs also state the invariants a
reader might otherwise have to discover by experiment (capacity vs length; "Vec
will never automatically shrink itself"; no small-value stack optimization). The
model to copy: contract + panics + cost + example, uniformly, per method. -->

```text
Vec<T> - a contiguous growable array type, written as Vec<T>, short for "vector".

Capacity and reallocation
  The capacity of a vector is the amount of space allocated for any future
  elements that will be added onto the vector. This is not to be confused with
  the length, which is the number of elements actually in the vector. If a
  vector's length exceeds its capacity, its capacity will automatically be
  increased, but its elements will have to be reallocated.

  For example, a vector with capacity 10 and length 0 would be an empty vector
  with space for 10 more elements. Pushing 10 or fewer elements will not change
  its capacity or cause reallocation. If the length is increased to 11, it will
  have to reallocate, which can be slow. For this reason it is recommended to
  use Vec::with_capacity whenever possible.

  Vec will never automatically shrink itself, even if completely empty.

pub fn push(&mut self, value: T)
  Appends an element to the back of a collection.

  # Panics
  Panics if the new capacity exceeds isize::MAX bytes.

  # Examples
  let mut vec = vec![1, 2];
  vec.push(3);
  assert_eq!(vec, [1, 2, 3]);

  # Time complexity
  Takes amortized O(1) time. If the vector's length would exceed its capacity
  after the push, O(capacity) time is taken to copy the vector's elements to a
  larger allocation.
```

---

## Negative pair (bad -> good)

<!-- Original, written for this corpus. The AI-reference failure mode is writing a
tutorial-flavored blurb that explains why you'd use the function instead of
stating its contract, and silently omitting error/edge behavior - exactly what
the two exemplars above refuse to do. -->

**BAD** - motivational prose, vague contract, no error cases:

```text
### fetchUser

This handy function lets you easily get a user from the database. Just pass
in an ID and it'll return the user for you. Great for when you need user data
in your app! Very fast and reliable.
```

**GOOD** - signature-first, exact contract, errors/panics/edges enumerated
(the SQLite + Rust discipline applied to an ordinary app function):

```text
### fetchUser(id: string, opts?: FetchOpts): Promise<User>

Fetches a single user by primary key.

Parameters:
  id                 Primary key. Must be a non-empty string; a v4 UUID by
                     convention but not validated as one.
  opts.timeoutMs     Per-call timeout. Default 5000. Must be > 0.
  opts.consistency   "strong" | "eventual". Default "strong".

Returns: a Promise resolving to the User record.

Throws / rejects:
  NotFoundError      no user with that id exists.
  ValidationError    id is empty or not a string.
  TimeoutError       no response within timeoutMs.
  ConnectionError    the pool could not acquire a connection.

Notes:
  - Does NOT cache. Each call hits the database.
  - Safe to call concurrently; uses a connection from the shared pool.
  - With consistency:"eventual" a read may return a user up to the replica
    lag (typically < 1s) out of date.

Example:
    const user = await fetchUser("a1b2-...", { timeoutMs: 2000 });

See also: fetchUsers(), User, FetchOpts.
```

<!-- Why the good one works: signature on the first line; every parameter has its
type, default, and constraint; every failure mode is a named error type rather
than "might fail"; the non-obvious behaviors (no caching, concurrency-safe,
staleness bound) are stated so the reader never has to test to find out; and the
example shows the one call, not a workflow. This is the sqlite3_step return-code
table and the Rust "# Panics/# Examples" discipline transferred to app code. -->