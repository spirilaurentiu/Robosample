<!--
CATEGORY: review (code review comments and exchanges - giving and receiving)

Exemplars in this file and their provenance:
  - Google eng-practices "How to Do a Code Review" · google/eng-practices
      · CC-BY 3.0 - excerpts with attribution (verify repo LICENSE before shipping)
  - Chromium "Respectful Code Reviews" (cr_respect.md) · chromium/src/docs
      · BSD-3-Clause (Chromium license) - excerpts with attribution
Link-only (individually authored, NO blanket license - never vendor):
  - LKML threads      https://lore.kernel.org/lkml/   (curate constructive
    threads only; the famous flame wars are anti-exemplars, not style guides)
  - GitHub PR comments, Gerrit review threads - same status: each comment is
    copyright its author, no license grant, regardless of the repo's license.
The worked exchange at the bottom is ORIGINAL FICTION for this corpus (clearly
marked) - kept because a full both-sides exchange with severity labels is hard
to source under license, but it is a supplement to the real texts, not their
substitute.
-->

# Code review

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. The standard is IMPROVEMENT, not perfection: approve once the change
   definitely improves overall code health, even if it isn't perfect. Don't
   block on preferences; do block on correctness and health regressions.
2. Review the right things in the right order: design first (does this change
   belong, does it fit the system?), then functionality (does it do what the
   author intended, and is that good for users?), then complexity, then tests -
   style last, and only where a style guide backs it.
3. Every blocking comment states a REASON grounded in principles, data, or the
   codebase - not "I wouldn't do it this way." If several approaches are validly
   equal, the author's choice stands.
4. Label severity so the author can triage: blocking issue vs. "Nit:" vs.
   "Optional/Consider:". Unlabeled comments all read as demands.
5. Comment on the code, never the coder ("this function re-hashes on every
   part" not "you got the hashing wrong"). Ask questions where you might be
   missing context; the author usually knows something you don't.
6. Say what's good, specifically. Reinforcing a pattern worth repeating is
   review signal too, not politeness filler.
7. As the AUTHOR: respond to every comment, answer the question not just the
   letter ("Done" vs. explaining), push back with reasons when you disagree,
   and never take it personally - review is of the code.
8. Watch for over-engineering: solving speculative future problems is a
   complexity cost now. Encourage solving the problem that exists.
Anti-patterns: "LGTM" on a change you didn't understand; a wall of unlabeled
nits burying one real defect; blocking on personal preference; sarcasm;
rewriting the patch in comments instead of stating the requirement; approving
to avoid conflict.
-->

## Google eng-practices - the reviewer's checklist (CC-BY, attributed)

<!-- Annotation: The canonical statement of WHAT to review and in what priority
order. Note the framing of functionality as two distinct questions - does the
code do what the author intended, AND is that behavior good for users - and the
explicit vigilance against over-engineering: reviewers are told to push back on
generality nobody needs yet. Attribution: Google eng-practices,
https://google.github.io/eng-practices/ , CC-BY. -->

```text
In doing a code review, you should make sure that:

  Design:        Do the interactions of various pieces of code in the CL make
                 sense? Does this change belong in your codebase, or in a
                 library? Does it integrate well with the rest of your system?
                 Is now a good time to add this functionality?
  Functionality: Does the code behave as the author likely intended? Is the
                 way the code behaves good for its users?
  Complexity:    Could the code be made simpler? Would another developer be
                 able to easily understand and use this code when they come
                 across it in the future?
  Tests:         Does the code have correct and well-designed automated tests?

A particular type of complexity is over-engineering, where developers have made
the code more generic than it needs to be, or added functionality that isn't
presently needed by the system. Reviewers should be especially vigilant about
over-engineering. Encourage developers to solve the problem they know needs to
be solved now, not the problem that the developer speculates might need to be
solved in the future.

On conflict: the first step should always be for the developer and reviewer to
try to come to consensus. Review decisions are based on underlying principles,
not personal opinion. If the author can demonstrate (through data or solid
engineering principles) that several approaches are equally valid, the reviewer
should accept the preference of the author. Don't let a CL sit around because
the author and the reviewer can't come to an agreement - escalate.
```

## Chromium - "Respectful Code Reviews" (BSD-3-Clause, attributed)

<!-- Annotation: Real text by Chromium engineers, and the best licensed source of
HOW TO SAY IT - it gives actual phrasings to avoid and their replacements, plus
the epistemics behind the tone: assume the disagreement comes from an
information gap, not incompetence, so a review comment's job is to transfer the
missing information. Attribution: Chromium project, docs/cr_respect.md,
https://chromium.googlesource.com/chromium/src/+/main/docs/cr_respect.md ,
BSD-3-Clause. -->

```text
We attract competent people - and that means even when they're wrong, it most
likely comes from lack of information, not from inability. A "bad" CL usually
means one of the parties is in possession of information the other one isn't
aware of.

It might be obvious to you that some code is wrong, but it's probably not
obvious to the author - or they wouldn't have written it that way. So please
don't say "This is wrong". Instead, explain at least what the right way looks
like.

Please don't say things like "no sane person would ever do this" or "this
algorithm is terrible", whether it's about the change you're reviewing or about
the surrounding code. While it might intimidate the reviewee into doing what
you want, it's not helpful in the long run - they will feel incapable, and
there is not much info in there to help them improve. "This is a good start,
but it could use some work" or "This needs some cleanup" are nicer ways of
saying it. Discuss the code, not the person.

If there is a disagreement, have a quick in-person/video/IM chat to sort out
what is going on - it's much easier to address all the little "Oh, I didn't
know"s in a single face-to-face than back-and-forth with long delays. And
please make sure to record the outcomes on the review.
```

## Worked review exchange (ORIGINAL FICTION for this corpus - supplement, not source)

<!-- Original, written for this corpus, and FICTIONAL - continues the multipart
upload bug: the author submits the `??=` fix that issues.md's root-cause trace
produced. Shows both roles: a reviewer who catches a real gap without blocking on
noise, labels severities, and reinforces what's good; an author who answers the
substance, pushes back once with a reason, and concedes where the reviewer is
right. -->

```text
PR: fix multipart uploads signing every part with part 1's payload hash
    (fixes #482 - SignatureDoesNotMatch on files > 5 MB)

  - ctx.headers["x-amz-content-sha256"] ??= sha256(part);
  + ctx.headers["x-amz-content-sha256"] = sha256(part);
  + test: two parts signed through one context get different payload hashes

REVIEWER:
  Nice find, and the regression test pinning two-parts-differ is exactly the
  right shape - it fails on the old code for the right reason.

  [blocking] The fix corrects the hash, but the underlying hazard remains: the
  reused context still lets ANY per-part field silently leak across parts -
  we're one more `??=` away from the same class of bug. Can we make the
  boundary structural instead? E.g. sign from `{...ctx, ...perPartFields}` so
  the reused object simply never holds per-part state. Grounds: this is the
  second cross-part leak this quarter (#461 was retry headers).

  [question] Does anything else currently piggyback on the context between
  parts - content-type, checksums? If yes they have the same latent bug and
  this PR is the moment to catch them.

  Nit: test name says "signs parts independently"; it specifically asserts
  payload hashes differ - worth saying that.

AUTHOR:
  Good point on the class-not-instance fix - done. Per-part fields now live in
  a fresh object merged at sign time; the shared ctx is frozen after creation,
  so the next `??=` is a TypeError in dev instead of a silent leak.

  Checked the other headers: content-type is constant across parts (fine);
  x-amz-checksum-crc32 had the SAME latent bug - never reported because we
  don't enable checksums by default. Fixed + test added. Good catch.

  On freezing in production too: pushing back - Object.freeze on the hot path
  costs ~3% in the upload benchmark and the dev-mode freeze already catches
  misuse before it ships. Keeping it dev-only unless you feel strongly.

  Nit: renamed.

REVIEWER:
  Dev-only freeze + the structural merge is the right trade, agreed - the 3%
  matters more than belt-and-suspenders. LGTM. Thanks for chasing the checksum
  variant; that one would have bitten us the day someone flips the flag.
```

<!-- Why this exchange works: the reviewer's blocking comment targets the CLASS
of defect, not the instance, and grounds it in evidence (a prior incident), so
it reads as engineering, not taste; severities are labeled so the author can
triage; the question surfaces a second latent bug the reviewer suspected but
couldn't see; the author answers substance rather than "Done"-ing everything,
pushes back exactly once WITH data (the benchmark), and the reviewer accepts -
modelling that authors win arguments they can support. Approval lands when the
change improves code health, not when it's perfect. -->

---

## Negative pair (bad -> good) - the same defect, two review comments

<!-- Original. The failure mode isn't rudeness alone - it's comments that give
the author nothing to act on: no location, no mechanism, no requirement, no
severity, coder-directed instead of code-directed. -->

**BAD** - vague, personal, no mechanism, no requirement:

```text
This caching stuff is wrong and honestly kind of sloppy. Did you even test
multipart? Rewrite the signing code properly.
```

**GOOD** - located, mechanism stated, requirement explicit, severity labeled, code-directed:

```text
[blocking] multipart.js:86 - `??=` computes the payload hash only when the
field is unset, but since a3f19c2 the context is reused across parts, so parts
2..n get signed with part 1's hash. That's the #482 SignatureDoesNotMatch.
Needs: per-part hash computed unconditionally, plus a test asserting two parts
signed through one context produce different x-amz-content-sha256 values -
otherwise this regresses the next time someone "optimizes" the header setup.
```

<!-- Why the good one works: file:line, the exact mechanism (what `??=` does
under context reuse), the observed consequence tied to the tracked issue, a
concrete acceptance requirement including the test that pins it, and a stated
severity - everything the author needs to fix it without a follow-up round.
And not one word about the author. -->