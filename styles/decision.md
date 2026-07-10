<!--
CATEGORY: decision (design proposals / rationale documents)

Exemplars in this file and their provenance:
  - Swift SE-0269 - swiftlang/swift-evolution - Apache-2.0 WITH LLVM-exception
  - Rust RFC 0169 - rust-lang/rfcs            - MIT OR Apache-2.0
-->

# Decision

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. Motivation is concrete pain, not assertion. Show the friction with a before
   snippet (SE-0269 opens on real code littered with redundant `self.`).
2. Walk the arc: Motivation -> Proposed solution -> Detailed design ->
   Alternatives considered -> (Future directions). Don't skip alternatives.
3. Reject alternatives *with reasons*, and record honest uncertainty
   (Rust RFC 169: "I cannot figure out if this is genius or madness").
4. Scope explicitly. Name what you are deliberately NOT doing and defer it
   (SE-0269 pushes `weak self` handling to Future Directions with rationale).
5. State impact plainly (source compatibility, ABI) even when the answer is "none".
Anti-patterns: a solution with no alternatives section; "this is obviously
better"; motivation that asserts importance instead of showing cost.
-->

## Swift proposal `SE-0269` - Increase availability of implicit `self` in `@escaping` closures

<!-- Annotation: A strong, ACCEPTED proposal (Implemented, Swift 5.3). Note the
Motivation: it doesn't claim the feature is important, it shows redundant `self.`
piling up in real async code and explains why the first explicit `self` already
satisfies the compiler's intent-signalling goal. Alternatives Considered lists
the more aggressive option (always require `self` in the capture list) and says
why it wasn't pursued (source-breaking). Future Directions fences off weak
captures deliberately rather than half-solving them. -->

````text
# Increase availability of implicit `self` in `@escaping` closures when reference cycles are unlikely to occur

* Proposal: [SE-0269](0269-implicit-self-explicit-capture.md)
* Author: [Frederick Kellison-Linn](https://github.com/jumhyn)
* Review Manager: [Ted Kremenek](https://github.com/tkremenek)
* Status: **Implemented (Swift 5.3)**
* Implementation: [apple/swift#23934](https://github.com/apple/swift/pull/23934)
* Bug: [SR-10218](https://bugs.swift.org/browse/SR-10218)
* Forum threads: [Discussion](https://forums.swift.org/t/review-capture-semantics-of-self/22017), [Pitch](https://forums.swift.org/t/allow-implicit-self-in-escaping-closures-when-self-is-explicitly-captured/22590), [Review thread](https://forums.swift.org/t/se-0269-increase-availability-of-implicit-self-in-escaping-closures-when-reference-cycles-are-unlikely-to-occur/30376)

## Introduction

Modify the rule that all uses of `self` in escaping closures must be explicit by allowing for implicit uses of `self` in situations where the user has already made their intent explicit, or where strong reference cycles are otherwise unlikely to occur. There are two situations covered by this proposal. The first is when the user has explicitly captured `self` in the closure's capture list, so that the following would compile without error:

```swift
class Test {
    var x = 0
    func execute(_ work: @escaping () -> Void) {
        work()
    }
    func method() {
        execute { [self] in
            x += 1
        }
    }
}
```

Secondly, this proposal would make implicit `self` available in escaping closures when `self` is a value type, so that the following would become valid:

```swift
struct Test {
    var x = 0
    func execute(_ work: @escaping () -> Void) {
        work()
    }
    func method() {
        execute {
            print(x)
        }
    }
}
```

## Motivation

In order to prevent users from inadvertently creating retain cycles, the Swift compiler today requires all uses of  `self` in escaping closures to be explicit. Attempting to reference a member `x` of `self` without the `self` keyword gives the error:

```
error: reference to property 'x' in closure requires explicit 'self.' to make capture semantics explicit
```

In codebases that choose to omit `self` where possible, this can result in a lot of unwanted noise, if many properties of `self` are used in a row:

```swift
execute {
    let foo = self.doFirstThing()
    performWork(with: self.bar)
    self.doSecondThing(with: foo)
    self.cleanup()
}
```

It also results in a lot of unnecessary repetition. The motivation for requiring explicit usage of `self` is to force the user to make the intent to capture `self` explicit, but that goal is accomplished after the first explicit usage of `self` and is not furthered by any of the subsequent usages.

In codebases that make heavy use of asynchronous code, such as clients of [PromiseKit](https://github.com/mxcl/PromiseKit), or even Apple's new [Combine](https://developer.apple.com/documentation/combine) and [SwiftUI](https://developer.apple.com/documentation/swiftui) libraries, maintainers must either choose to adopt style guides which require the use of explicit `self` in all cases, or else resign to the reality that explicit usage of `self` will appear inconsistently throughout the code. Given that Apple's preferred/recommended style is to omit `self` where possible (as evidenced by examples throughout [The Swift Programming Language](https://docs.swift.org/swift-book/) and the [SwiftUI Tutorials](https://developer.apple.com/tutorials/swiftui)), having any asynchronous code littered with `self.` is a suboptimal state of affairs.

The error is also overly conservative-it will prevent the usage of implicit `self` even when it is very unlikely that the capture of `self` will somehow cause a reference cycle. With function builders in SwiftUI, the automatic resolution to the "make capture semantics explicit" error (and indeed, the fix-it that is currently supplied to the programmer) is just to slap a `self.` on wherever the compiler tells you to. While this is likely fine when `self` is a value type, building this habit could cause users to ignore the error and apply the fix-it in cases where it is not in fact appropriate.

There are also cases that are just as (if not more) likely to cause reference cycles that the tool currently misses. Once such instance is discussed in **Future Directions** below. These false negatives are not addressed in this proposal, but improving the precision of this diagnostic will make the tool more powerful and less likely to be dismissed if (or when) new diagnostic cases are introduced.

## Proposed solution

First, allow the use of implicit `self` when it appears in the closure's capture list. The above code could then be written as:

```swift
execute { [self] in
    let foo = doFirstThing()
    performWork(with: bar)
    doSecondThing(with: foo)
    cleanup()
}
```

This change still forces code which captures `self` to be explicit about its intentions, but reduces the cost of that explicitness to a single declaration. With this change explicit capture of `self` would be one of *two* ways to get rid of the error, with the current method of adding `self.` to each property/method access (without adding `self` to the capture list) remaining as a viable option.

The compiler would also offer an additional fix-it when implicit `self` is used:

```swift
execute { // <- Fix-it: capture 'self' explicitly to enable implicit 'self' in this closure. Fix-it: insert '[self] in'
    let foo = doFirstThing()
    performWork(with: bar)
    doSecondThing(with: foo)
    cleanup()
}
```

Second, if `self` is a value type, we will not require any explicit usage of `self` (at the call/use site or in the capture list), so that if `self` were a `struct` or `enum` then the above could be written as simply:

```swift
execute {
    let foo = doFirstThing()
    performWork(with: bar)
    doSecondThing(with: foo)
    cleanup()
}
```

## Detailed design

Whenever `self` is declared explicitly in an escaping closure's capture list, or its type is a value type, any code inside that closure can use names which resolve to members of the enclosing type, without specifying `self.` explicitly. In nested closures, the *innermost* escaping closure must capture `self`, so the following code would be invalid:

```swift
execute { [self] in
    execute { // Fix-it: capture 'self' explicitly to enable implicit 'self' in this closure.
        x += 1 // Error: reference to property 'x' in closure requires explicit use of 'self' to make capture semantics explicit. Fix-it: reference 'self.' explicitly.
    }
}
```

This new behavior will also be available for closures with  `unowned(safe)` and `unowned(unsafe)` captures of `self` since the programmer has clearly declared their intent to capture `self` (as for strong captures). For `weak` captures of `self`, it's not immediately clear what bare reference to an instance property even *means* when the local `self` is bound weakly-should it be treated as though the programmer had written `self?.`? Should there be a special syntax available in this context, such as `?.x`? Should we offer a diagnostic offering to insert `self?.`? Or should we instead suggest that the programmer re-bind `self` strongly via a `guard` statement? There are enough open questions here that have not been sufficiently discussed, so this proposal leaves the handling of `weak self` captures as a future direction worthy of further consideration.

The existing errors and fix-its have their language updated accordingly to indicate that there are now multiple ways to resolve the error. In addition to the changes noted above, we will also have:

```
Error: call to method <method name> in closure requires explicit use of 'self' to make capture semantics explicit.
```

The new fix-it for explicitly capturing `self` will take slightly different forms depending on whether there is a capture list already present ("insert '`self, `'"), whether there are explicit parameters ("insert '`[self]`'"), and whether the user has already captured a variable with the name `self` (in which case the fix-it would not be offered). Since empty capture lists are allowed (`{ [] in ... }`), the fix-it will cover this case too.

This new rule would only apply when `self` is captured directly, and with the name `self`. This includes captures of the form `[self = self]` but would still not permit implicit `self` if the capture were `[y = self]`. In the unusual event that the user has captured another variable with the name `self` (e.g. `[self = "hello"]`), we will offer a note that this does not enable use of implicit `self` (in addition to the existing error attached to the attempted use of implicit `self`):

```swift
Note: variable other than 'self' captured here under the name 'self' does not enable implicit 'self'.
```

If the user has a capture of `weak self` already, we offer a special diagnostic note indicating that `weak` captures of `self` do not enable implicit `self`:

```swift
Note: weak capture of 'self' here does not enable implicit 'self'.
```

If either of the two above notes is present, we will not offer the usual fix-its for resolving this error, since the code inserted would be erroneous.

## Source compatibility

This proposal makes previously illegal syntax legal, and has no effect on source compatibility.

## Effect on ABI stability

This is an entirely frontend change and has no effect on ABI stability.

## Effect on API resilience

This is not an API-level change, and has no effect on API resilience.

## Future Directions

### Bound method references

While this proposal opens up implicit `self` in situations where we can be reasonably sure that we will not cause a reference cycle, there are other cases where implicit `self` is currently allowed that we may want to disallow in the future. One of these is allowing bound method references to be passed into escaping contexts without making it clear that such a reference captures `self`. For example:

```swift
class Test {
    var x = 0
    func execute(_ work: @escaping () -> Void) {
        work()
    }
    func method() {
        execute(inc) // self is captured, but no error!
        execute { inc() } // error
    }
    func inc() {
        x += 1
    }
}
```

This would help reduce the false negatives of the "make capture semantics explicit" error, making it more useful and hopefully catching reference cycles that the programmer was previously unaware of.

### Weak captures of `self`

It might be desirable to add improved diagnostic information (or similarly do away with the error) in the case where the programmer has explicitly captured `self`, but only as a weak reference. For the reasons noted in this proposal, this was not pursued, but answering the lingering questions and issuing a follow-up proposal would fit nicely with this direction.

## Alternatives considered

### Always require `self` in the capture list

The rule requiring the use of explicit `self` is helpful when the code base does not already require `self` to be used on all instance accesses. However, many code bases have linters or style guides that require `self` to be used explicitly always, making the capture semantics opaque. Always requiring `self` to be captured in the capture list explicitly would ensure that there are no `self` captures that the programmer is unaware of, even if they naturally use `self` for instance accesses. This would be a more drastic, source breaking change (and is not ruled out by adopting this change), so it was not seriously pursued as part of this proposal.

### Eliminate the former fix-it

A less extreme solution to the problem described above is to simply stop offering the current fix-it that suggests adding the explicit `self.` at the point of reference in favor of only recommending the explicit capture list fix-it, when possible.
````

## Rust RFC `0169` - `use PATH as ID` rebinding syntax

<!-- Annotation: The whole decision arc in miniature - a small, well-argued change.
Motivation is visual and concrete (imports don't line up; your eye scans both
sides). Alternatives enumerates every rejected token (@, =, :) with a one-line
reason each, including self-aware doubt. Ends with Unresolved Questions rather
than pretending everything is settled. If you keep one decision exemplar, keep
this one. -->

```text

    Start Date: 2014-07-16
    RFC PR #: #169
    Rust Issue #: https://github.com/rust-lang/rust/issues/16461

Summary

Change the rebinding syntax from use ID = PATH to use PATH as ID, so that paths all line up on the left side, and imported identifiers are all on the right side. Also modify extern crate syntax analogously, for consistency.
Motivation

Currently, the view items at the start of a module look something like this:

mod old_code {
  use a::b::c::d::www;
  use a::b::c::e::xxx;
  use yyy = a::b::yummy;
  use a::b::c::g::zzz;
}

This means that if you want to see what identifiers have been imported, your eyes need to scan back and forth on both the left-hand side (immediately beside the use) and the right-hand side (at the end of each line). In particular, note that yummy is not in scope within the body of old_code

This RFC proposes changing the grammar of Rust so that the example above would look like this:

mod new_code {
  use a::b::c::d::www;
  use a::b::c::e::xxx;
  use a::b::yummy as yyy;
  use a::b::c::g::zzz;
}

There are two benefits we can see by comparing mod old_code and mod new_code:

    As alluded to above, now all of the imported identifiers are on the right-hand side of the block of view items.

    Additionally, the left-hand side looks much more regular, since one sees the straight lines of a::b:: characters all the way down, which makes the actual differences between the different paths more visually apparent.

Detailed design

Currently, the grammar for use statements is something like:

  use_decl : "pub" ? "use" [ ident '=' path
                            | path_glob ] ;

Likewise, the grammar for extern crate declarations is something like:

  extern_crate_decl : "extern" "crate" ident [ '(' link_attrs ')' ] ? [ '=' string_lit ] ? ;

This RFC proposes changing the grammar for use statements to something like:

  use_decl : "pub" ? "use" [ path "as" ident
                            | path_glob ] ;

and the grammar for extern crate declarations to something like:

  extern_crate_decl : "extern" "crate" [ string_lit "as" ] ? ident [ '(' link_attrs ')' ] ? ;

Both use and pub use forms are changed to use path as ident instead of ident = path. The form use path as ident has the same constraints and meaning that use ident = path has today.

Nothing about path globs is changed; the view items that use ident = path are disjoint from the view items that use path globs, and that continues to be the case under path as ident.

The old syntaxes "use" ident '=' path and "extern" "crate" ident '=' string_lit are removed (or at least deprecated).
Drawbacks

    pub use export = import_path may be preferred over pub use import_path as export since people are used to seeing the name exported by a pub item on the left-hand side of an = sign. (See “Have distinct rebinding syntaxes for use and pub use” below.)

    The ‘as’ keyword is not currently used for any binding form in Rust. Adopting this RFC would change that precedent. (See “Change the signaling token” below.)

Alternatives
Keep things as they are

This just has the drawbacks outlined in the motivation: the left-hand side of the view items are less regular, and one needs to scan both the left- and right-hand sides to see all the imported identifiers.
Change the signaling token

Go ahead with switch, so imported identifier is on the left-hand side, but use a different token than as to signal a rebinding.

For example, we could use @, as an analogy with its use as a binding operator in match expressions:

mod new_code {
  use a::b::c::d::www;
  use a::b::c::e::xxx;
  use a::b::yummy @ yyy;
  use a::b::c::g::zzz;
}

(I do not object to path @ ident, though I find it somehow more “line-noisy” than as in this context.)

Or, we could use =:

mod new_code {
  use a::b::c::d::www;
  use a::b::c::e::xxx;
  use a::b::yummy = yyy;
  use a::b::c::g::zzz;
}

(I do object to path = ident, since typically when = is used to bind, the identifier being bound occurs on the left-hand side.)

Or, we could use :, by (weak) analogy with struct pattern syntax:

mod new_code {
  use a::b::c::d::www;
  use a::b::c::e::xxx;
  use a::b::yummy : yyy;
  use a::b::c::g::zzz;
}

(I cannot figure out if this is genius or madness. Probably madness, especially if one is allowed to omit the whitespace around the :)
Have distinct rebinding syntaxes for use and pub use

If people really like having ident = path for pub use, by the reasoning presented above that people are used to seeing the name exported by a pub item on the left-hand side of an = sign, then we could support that by continuing to support pub use ident = path.

If we were to go down that route, I would prefer to have distinct notions of the exported name and imported name, so that:

pub use a = foo::bar; would actually import bar (and a would just be visible as an export), and then one could rebind for export and import simultaneously, like so: pub use exported_bar = foo::bar as imported_bar;

But really, is pub use foo::bar as a all that bad?
Allow extern crate ident as ident

As written, this RFC allows for two variants of extern_crate_decl:

extern crate old_name;
extern crate "old_name" as new_name;

These are just analogous to the current options that use = instead of as.

However, the RFC comment dialogue suggested also allowing a renaming form that does not use a string literal:

extern crate old_name as new_name;

I have no opinion on whether this should be added or not. Arguably this choice is orthogonal to the goals of this RFC (since, if this is a good idea, it could just as well be implemented with the = syntax). Perhaps it should just be filed as a separate RFC on its own.
Unresolved questions

    In the revised extern crate form, is it best to put the link_attrs after the identifier, as written above? Or would it be better for them to come after the string_literal when using the extern crate string_literal as ident form?
```
