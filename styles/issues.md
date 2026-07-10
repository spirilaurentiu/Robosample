<!--
CATEGORY: issues (bug reports, feature requests, and triage-ready problem writeups)

Exemplars in this file and their provenance:
  - Rubric + both worked examples below - original for this corpus - no license

A good bug report IS a completed investigation: to write "reproduces 5/5,
minimal repro is X, expected Y got Z, hypothesis fenced off," you must have
actually reproduced, minimized, and localized. So the report rubric doubles
as an INVESTIGATION CHECKLIST for the reproduce -> minimize -> localize loop.
The "Investigation trace" section below teaches that loop explicitly.

This category covers the FULL loop, in two phases with different rules:
  Phase 1 (reporting):    reproduce -> minimize -> localize to a code path ->
    labelled hypothesis. Output: the report. Theories stay fenced off.
  Phase 2 (root-causing): hypothesis -> falsifying experiments -> bisect ->
    the buggy line -> fix verified against the phase-1 repro. Output: the
    root-cause trace. Here theories are the tool - but each one must be
    TESTED before it's believed, and killed when the evidence says so.
The phase-1 report is the contract between the two: its minimal repro is the
test harness for every phase-2 experiment and the final proof of the fix.
-->

# Issues / bug reports

<!-- RUBRIC - distilled from Simon Tatham's essay and curl's guide; what a good report contains
1. The goal is REPRODUCIBILITY: give the reader everything they need to make the
   bug happen on their own machine. That is the whole job.
2. Separate three things and keep them separate:
     - what you did (exact steps / commands, verbatim),
     - what you expected,
     - what actually happened (verbatim output, not a paraphrase).
3. Report symptoms, not your theory of the cause. A diagnosis stated as fact
   sends the maintainer down your wrong path. If you have a hypothesis, label it.
4. Include environment: exact versions of the program, the OS, and anything else
   in the loop.
5. Minimize: reduce to the smallest input/steps that still trigger it. A minimal
   repro is the single most valuable thing you can provide.
6. State reproducibility rate honestly ("5/5" vs "intermittent, ~3/10").
Anti-patterns: "it doesn't work"; pasting a fix opinion instead of the symptom;
screenshots of text that should be copy-pasted; no version, no repro steps.
-->

## Investigation trace (the bug-hunting loop, then the report it produces)

<!-- Original, written for this corpus. This teaches the HUNT, not just the writeup:
the same discipline that makes a good report (reproduce, minimize, localize, form
a labelled hypothesis) is what finds the bug. An agent can be handed a symptom +
codebase and asked to produce this trace; the final report falls out of it. -->

```text
SYMPTOM (as first observed)
  "Uploads sometimes fail." Vague. Not yet reproducible, not yet minimal.

STEP 1 - Reproduce reliably.
  Try the obvious path. A 2 MB file: succeeds, 5/5. A 50 MB file: fails, 5/5.
  => Not random. Size-dependent. Reproducibility went from "sometimes" to "5/5
     above a threshold" - the single most important step.

STEP 2 - Bisect the input to find the boundary.
  Binary-search the file size:
     4 MB ok, 8 MB fail; 5 MB ok, 6 MB fail; narrow to ~5 MB.
  => Threshold is exactly 5 MB - which is this SDK's multipart-upload cutoff.
     The single-part path works; the multipart path fails. Search space just
     collapsed from "the uploader" to "the multipart path."

STEP 3 - Minimize the repro.
  Smallest trigger: one uploadFile() call on any file > 5 MB. Drop everything
  else (no retries, no concurrency, fresh client). Still fails 5/5.
  => Minimal repro is 3 lines. Rules out retry/race theories entirely.

STEP 4 - Localize (as far as a reporter should go).
  The error is "SignatureDoesNotMatch", and it only appears on the multipart
  path. That points at part-signing, not transport. That is enough to file:
  a labelled hypothesis, not a claimed root cause.

STEP 5 - Write it up. The report below is just this trace, cleaned up.
```

## The report that trace produces (bad → good)

**BAD** - theory-as-fact, no repro, no environment (what the symptom looked like before the hunt):

```text
Title: Upload is broken

The S3 upload is completely broken, I think there's a race condition in the
retry logic. It keeps failing for me. Can you fix the mutex?
```

**GOOD** - symptom-first, reproducible, minimized, environment stated, theory labelled:

```text
Title: uploadFile() throws "SignatureDoesNotMatch" on files > 5 MB (multipart path)

Environment:
  - foo-uploader 2.4.1
  - Node 20.11.0, macOS 14.5 (arm64)
  - AWS SDK v3.577.0, region eu-central-1

What I did (reproduces 5/5):
  1. const up = new Uploader({ region: "eu-central-1" });
  2. await up.uploadFile("./10mb.bin", "my-bucket", "test.bin");
     (any file over the 5 MB multipart threshold)

Expected: file uploaded, resolves with the object URL.

Actual: rejects after ~2s with:
    S3ServiceException: SignatureDoesNotMatch: The request signature we
    calculated does not match the signature you provided.
    at ... uploader/multipart.js:88

Scope: files <= 5 MB (single-part path) upload fine, 5/5. Only the multipart
path (> 5 MB) fails. Minimal repro attached (10mb.bin from
`head -c 10m /dev/urandom`).

Hypothesis (unconfirmed): part-signing may hash the buffer after the stream has
consumed it, but I haven't verified - the symptom above is what's certain.
```

<!-- Why the good one works: the title is a specific, searchable symptom with the
exact threshold; environment has exact versions; the steps are copy-pasteable and
state a reproduction rate; the actual output is verbatim including error class and
file:line; the single-part-vs-multipart scope line already narrows the search to
one code path (the payoff of STEP 2); and the theory is present but explicitly
fenced off as unconfirmed so it informs without misdirecting. The report is the
investigation, written down. -->

---

## Phase 2 - Root-cause trace (from labelled hypothesis to the buggy line)

<!-- RUBRIC for root-causing - how phase 2 differs from phase 1
1. Every hypothesis gets an EXPERIMENT that could falsify it, run against the
   phase-1 minimal repro. No experiment, no belief.
2. Prefer experiments that split the search space, not ones that confirm your
   favorite theory: differential tests (works-here/fails-there), git bisect,
   instrumenting the boundary between two suspects.
3. Kill hypotheses out loud. A dead theory is progress - record it and why it
   died, so nobody (including future-you) re-walks that path.
4. "Found the line" is not the end. Three more steps: explain WHY that line
   produces exactly the observed symptom (mechanism), write a minimal failing
   test that pins it, and verify the fix turns the phase-1 repro green.
5. Ask why the defect survived until now (missing test? untested path? recent
   regression?) - one sentence, feeds regression prevention.
Anti-patterns: fixing the symptom where it appears instead of where it's caused;
declaring root cause from the first confirming observation; "fixed it, tests
pass" without running the original repro; stacking a second hypothesis on an
untested first one.
-->

<!-- Original, written for this corpus. Continues the SAME bug: the phase-1
report ended with "Hypothesis (unconfirmed): part-signing may hash the buffer
after the stream has consumed it." Phase 2 starts by trying to kill that
hypothesis - and follows the evidence somewhere else first. -->

```text
INPUT: the phase-1 report. Repro: uploadFile() on any file > 5 MB,
SignatureDoesNotMatch at multipart.js:88, 5/5. Hypothesis H1 (unconfirmed):
part-signing hashes the buffer after the stream consumed it.

STEP 1 - Turn the hypothesis into a falsifiable experiment.
  If H1 is true, signing should see an EMPTY/partial buffer. Test: log the
  byte length the signer hashes vs. the part length actually sent.
    signer hashed: 5242880 bytes   |   part sent: 5242880 bytes
  Lengths match on every part. H1 as stated is DEAD. Record it, move on.
  (Note what killing H1 bought us: the buffer is intact, so the mismatch must
  be in WHAT is hashed or WHICH signature is attached - smaller search space.)

STEP 2 - Differential experiment to split the remaining space.
  Upload a 12 MB file (3 parts). Which part fails?
    part 1: 200 OK   part 2: 403 SignatureDoesNotMatch
  First part always succeeds; failure starts at part 2. The defect involves
  per-part STATE, not signing in general. New hypothesis H2: something from
  part 1's request is being reused for part 2 (stale signature, stale headers,
  stale hash).

STEP 3 - Bisect history to catch the regression (independent evidence line).
  The reporter said 2.4.0 worked. git bisect with the phase-1 repro as the
  test script:
    $ git bisect start v2.4.1 v2.4.0
    $ git bisect run node repro.js
    ...
    a3f19c2 is the first bad commit
      "perf: reuse request context across parts to cut allocations"
  A perf commit that REUSES a per-request object across parts - exactly the
  shape H2 predicts. Two independent lines of evidence now point at the same
  place. That convergence, not either line alone, is what justifies belief.

STEP 4 - Read the code at the intersection; name the mechanism.
  multipart.js:82-90 after a3f19c2:
    82  const ctx = this._ctx ??= makeSigningContext(creds);   // reused
    ...
    86  ctx.headers["x-amz-content-sha256"] ??= sha256(part);  // <-- BUG
    88  const sig = signV4(ctx);
  Line 86 uses `??=`: it computes the payload hash only if unset. On part 1
  the field is empty and gets part 1's hash. The context is now reused (line
  82, the perf change), so on part 2 the field is ALREADY set - part 2 is
  signed with part 1's payload hash. Mechanism explains every observation:
  part 1 succeeds, every later part 403s, single-part path (fresh context per
  request) unaffected, threshold exactly at the multipart cutoff.

STEP 5 - Pin it with a minimal failing test, then fix at the cause.
  Test (fails on HEAD): sign two different parts through one context; assert
  their x-amz-content-sha256 values differ.
  Fix: hash is per-part state, not context state - compute it unconditionally
  per part (`=`, not `??=`), or exclude it from the reused context.

STEP 6 - Verify against the ORIGINAL repro, not just the new test.
  phase-1 repro: 10mb.bin uploads, 5/5. 12 MB / 3-part case: 5/5.
  New unit test green; kept as the regression guard.

STEP 7 - One sentence on why it escaped: the perf change had benchmarks but no
  multi-part correctness test; the suite's largest fixture was 1 MB - below the
  multipart threshold, so the reused-context path was never exercised.
```

<!-- Why this trace works as an exemplar: the first move is trying to FALSIFY the
inherited hypothesis, and it dies - modelling that theories are inputs to
experiments, not conclusions; each subsequent step is chosen to split the search
space (which part fails; which commit introduced it) rather than to confirm a
hunch; belief arrives only when two independent evidence lines converge; "found
the line" is followed by mechanism (why ??= + context reuse produces exactly
this symptom table), a pinning test, and verification against the ORIGINAL
phase-1 repro; and it closes with the why-it-escaped sentence that turns a bug
fix into a prevented class of bugs. -->