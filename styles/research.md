<!--
CATEGORY: research (empirical investigation writeups - the researcher's voice,
in three registers: DISCOVERY - formulating the question and the recall-maximizing
search vocabulary that finds the literature; REPORTING - enabling replication of
what you measured; and EXPOSITION - teaching an insight you earned. An autonomous
deep research agent needs all three: discovery to search, reporting for findings,
exposition for synthesis.)

Exemplars in this file and their provenance - ALL vendorable, this file is
complete and self-contained at vendor time:
  - "On the Relative Motion of the Earth and the Luminiferous Ether"
      Michelson & Morley, American Journal of Science, 1887 - Public domain
  - "Why Momentum Really Works" - Gabriel Goh, Distill, 2017 - CC-BY 4.0
      (Distill's standard license). Attribute as: Goh, "Why Momentum Really
      Works", Distill, 2017. doi:10.23915/distill.00006
  - spell.py - Peter Norvig, from the pytudes repository - MIT
      (repository LICENSE is MIT, and the file itself carries the grant:
      "Copyright (c) 2007-2016 Peter Norvig, MIT license"). Code and the
      measured output of running it are vendored; the accompanying essay's
      prose is NOT vendored - its lessons are carried in the original
      annotations below, in this corpus's own words.
  - "The Chemical History of a Candle" - Michael Faraday, 1861 - Public
      domain (Project Gutenberg #14474)

Relationship to siblings: issues.md is the DEBUGGING loop (a defect exists;
find it). research.md is the DISCOVERY loop (a question exists; answer it
with experiments). Same epistemics - hypotheses are inputs to experiments,
belief requires convergent evidence, dead theories are recorded - but a
different deliverable: new knowledge, quantified, with its limits stated.
-->

# Research writeups

<!-- RUBRIC - what an agent should reproduce from these exemplars
1. State the QUESTION and the expected outcome BEFORE the result. Michelson &
   Morley derive and publish the predicted displacement (0.4 fringe) from
   theory first; only then do they report what they saw. Norvig's spell.py
   docstring-level goal is the same shape: scoped targets, declared up front,
   so the reader can check the result against the pre-registered bar.
2. Report the miss without spin - the null result IS the result. M&M's
   observed displacement was under 1/20 of prediction, and they quantify
   exactly what that kills ("quite small enough entirely to refute Fresnel's
   explanation"). Norvig's measured 75%/68% lands under his 80-90% target
   and the number is printed anyway. The trust of the whole genre is built
   on writers who publish the numbers that embarrass their hopes.
3. Turn "we found nothing" into a QUANTIFIED upper bound. M&M don't say the
   ether wind wasn't detected; they say the relative velocity "is probably
   less than one sixth the earth's orbital velocity, and certainly less than
   one-fourth" - a null with error bars is knowledge, a shrug is not.
4. Meet the reader's existing mental model head-on and say precisely where it
   fails. Distill's opening: state the popular story (momentum = heavy ball
   rolling downhill), grant what it gets right, then "this standard story
   isn't wrong, but it fails to explain many important behaviors." The
   researcher-voice twin of SQLite's "a feature, not a bug" move in
   documentation.md.
5. Choose the SMALLEST system rich enough to answer the question, and justify
   the choice explicitly. Goh picks the convex quadratic because it is "rich
   enough to reproduce momentum's local dynamics in real problems, and yet
   simple enough to be understood in closed form." Norvig builds a toy and
   scopes it as one. Faraday picks a candle. A minimal-sufficient model is a
   claim of understanding; a kitchen-sink model is an evasion of it.
6. Evaluation hygiene is stated, not assumed: a development set you may look
   at, a held-out set you may not touch until the method is frozen, named as
   such in the writeup. Norvig's test harness (vendored below) encodes this
   as two separate spelltest calls; the fictional GOOD example below shows
   the prose form.
7. Error analysis works from VERBATIM failures, categorized by mechanism, and
   future work is derived from those observed failure classes - never from
   speculation about what might be nice to add. spell.py's verbose mode
   exists precisely to print the failure table that drives the next
   iteration.
8. State limits honestly, including the limits of the limits. Distill proves
   momentum cannot be beaten within its algorithm class, then immediately
   bounds the claim itself: the result "must not be taken literally, but
   spiritually" - it doesn't preclude reformulating the problem. M&M close by
   scoping their own null: "only the orbital motion of the earth is
   considered." Knowing what a negative result does NOT say is as
   load-bearing as the result.
9. Keep measured, derived, and conjectured claims visibly separate. This is
   issues.md's "hypothesis, labelled, fenced off" discipline promoted to a
   whole document: every sentence should be traceable to an experiment, a
   derivation, or an explicit "probably / it seems" marker - M&M's certainly-
   less-than vs probably-less-than is exactly this distinction, made twice
   in one sentence.
Anti-patterns: hedge fog ("the results may potentially suggest") - uncertainty
as verbal tic instead of a number; "significantly improves" with no baseline,
no magnitude, no definition of significant; notation before motivation;
silently tuning on the held-out set; an abstract that claims more than the
experiments support; "it is trivial/obvious that"; a null result reported as
failure-with-a-shrug instead of as an upper bound; declaring victory from
confirming runs without an experiment that could have falsified.
-->

## Discovery register - formulating the question and the search vocabulary

<!-- The query and the reformulation are an intermediate representation: they run
before any writeup exists, they steer every retrieval, and a human approves them
at Checkpoint A. Bias enters the pipeline here. This register governs how you
FORMULATE, not how you report, and is the single source for the query-vocabulary
rules; the scout and researcher agents reference it rather than restating it. -->

Retrieval vocabulary is the opposite discipline from expression vocabulary. A
writeup uses one precise term per concept; a search maximizes recall. The two
SHALL NOT be conflated - narrowing a query to your preferred term is how an entire
research programme that names the idea differently goes unseen.

1. Recall first. Cast wide and shallow before deep. Missing a renamed programme
   is the expensive error; a redundant hit is free.
2. Each synonym is another query. For every concept, search its canonical name,
   historical names, competing names, abbreviations, and the broader / narrower /
   neighboring terms. A newly found synonym is fed back as a new query.
3. Search equivalent formulations, not identical terminology. Translate the
   problem into neighboring disciplines (statistics, optimization, robotics,
   numerical analysis) and search their words for it.
4. A query MUST NOT encode its own answer. "does X improve Y" retrieves only
   confirmation. Every essential claim SHALL also be searched in the vocabulary
   of its opponents - competing terms, "failure of X", "X considered harmful".
   This is the query-side twin of the verifier's rule that a check MUST NOT
   encode the claim it tests.
5. Motivation before notation. Formulate the question from the phenomenon, not
   from a symbol you already picked.
6. Freeze the question, not the words. Checkpoint A freezes the problem statement.
   The terminology set stays append-only: new synonyms and competing terms are
   added throughout, and that growth is recall work, not drift.

Stop when a full pass adds no new terminology, no new canonical author, and no new
conceptual reading - a quantified saturation, not a feeling.

**BAD** - anchors on the user's words, presupposes its answer

    User asked: "does the torsional HMC move mix faster with a smaller timestep?"
    Query 1: torsional HMC timestep faster mixing
    Query 2: torsional HMC small timestep convergence
    Query 3: Robosample torsion move step size
    -> three phrasings of one term set, all presupposing "smaller is faster";
       misses "step-size dependence of HMC acceptance", "integrator stability
       limit", "optimal acceptance rate 0.65", and the robotics/optimization
       literature entirely. Low recall, confirmation-biased.

**GOOD** - phenomenon-first, synonym-expanded, includes the falsifying vocabulary

    Question (frozen): what sets the step size that maximizes sampling efficiency
    of the torsional HMC move, and does reducing it help or hurt?
    Concept map, each line a query family:
      HMC step size        <- leapfrog step, integrator timestep, dt
      acceptance vs dt      <- optimal acceptance rate, 0.651, tuning dt
      efficiency metric     <- ESS/gradient, autocorrelation time, mixing rate
      failure mode          <- integrator instability, energy drift, step-size limit
      neighboring fields    <- MALA/ULA step size (stats), trust-region step (optim)
    Opponent queries        <- "smaller timestep worse HMC", "too small step size
                               diffusive", cost-per-effective-sample vs dt
    -> wide recall, crosses fields, and searches the vocabulary that could refute
       "smaller is faster", not only the vocabulary that confirms it.

## Michelson & Morley, 1887 (public domain) - prediction first, honest null, quantified bound

<!-- Annotation: The reporting register's founding document. Three moves in
three sentences. (1) The prediction is computed and stated BEFORE the
observation - 0.4 fringe, derived from the theory under test, so the
experiment is falsifiable in print. (2) The observed miss is reported as a
RATIO to prediction (< 1/20, probably < 1/40) - and notice the two-tier
uncertainty language, "certainly less than" vs "probably less than": measured
bound and estimated bound, kept separate inside a single sentence. (3) The
null is converted through the theory (displacement ~ v^2) into a quantified
upper bound on the physical quantity, then into the specific theoretical
casualty (Fresnel's aberration explanation). Finally, the scope limit on
their own conclusion: only orbital motion was considered - the null does not
yet cover the solar system's own motion, and they say so and propose the
follow-up observation schedule. -->

```text
The distance D was about eleven meters, or 2x10^7 wave-lengths of yellow
light; hence the displacement to be expected was 0.4 fringe. The actual
displacement was certainly less than the twentieth part of this, and
probably less than the fortieth part. But since the displacement is
proportional to the square of the velocity, the relative velocity of the
earth and the ether is probably less than one sixth the earth's orbital
velocity, and certainly less than one-fourth.

In what precedes, only the orbital motion of the earth is considered.

It appears, from all that precedes, reasonably certain that if there be any
relative motion between the earth and the luminiferous ether, it must be
small; quite small enough entirely to refute Fresnel's explanation of
aberration.
```

<!-- Note also what the full paper spends most of its pages on: the
apparatus (the sandstone slab floating on mercury, the multiple-reflection
light path that bought them the factor-of-ten sensitivity their 1881 attempt
lacked) and the sources of error they engineered away. Methods detail in
proportion to how much the conclusion leans on it - a null result is only as
strong as the demonstrated sensitivity of the instrument that failed to see
the effect. -->

## Norvig - spell.py and its measured results (MIT) - the artifact, its harness, its numbers

<!-- Annotation: A complete empirical study small enough to vendor whole: the
method (the code), the pre-declared evaluation protocol (two test sets, one
for development, one final - encoded as two separate spelltest calls on two
separate files), and the instrument for error analysis (verbose mode prints
every failure with the word counts that explain WHY the wrong candidate won).
The measured results of running it: 75% of 270 correct at ~41 words/second on
the development set, 68% of 400 at ~35 wps on the final set - published
although they fall short of the author's stated 80-90% accuracy goal (speed
and brevity goals were met). The lesson is structural, not just moral: the
harness makes honesty mechanical. When the final set is a file you run once,
there is nothing to fudge. -->

```python
"""Spelling Corrector in Python 3.
Copyright (c) 2007-2016 Peter Norvig
MIT license: www.opensource.org/licenses/mit-license.php"""

import re
from collections import Counter

def words(text): return re.findall(r'\w+', text.lower())

WORDS = Counter(words(open('big.txt').read()))

def P(word, N=sum(WORDS.values())):
    "Probability of `word`."
    return WORDS[word] / N

def correction(word):
    "Most probable spelling correction for word."
    return max(candidates(word), key=P)

def candidates(word):
    "Generate possible spelling corrections for word."
    return (known([word]) or known(edits1(word)) or known(edits2(word)) or [word])

def known(words):
    "The subset of `words` that appear in the dictionary of WORDS."
    return set(w for w in words if w in WORDS)

def edits1(word):
    "All edits that are one edit away from `word`."
    letters    = 'abcdefghijklmnopqrstuvwxyz'
    splits     = [(word[:i], word[i:])    for i in range(len(word) + 1)]
    deletes    = [L + R[1:]               for L, R in splits if R]
    transposes = [L + R[1] + R[0] + R[2:] for L, R in splits if len(R)>1]
    replaces   = [L + c + R[1:]           for L, R in splits if R for c in letters]
    inserts    = [L + c + R               for L, R in splits for c in letters]
    return set(deletes + transposes + replaces + inserts)

def edits2(word):
    "All edits that are two edits away from `word`."
    return (e2 for e1 in edits1(word) for e2 in edits1(e1))

################ Test harness: dev set and a separate, final test set

def spelltest(tests, verbose=False):
    "Run correction(wrong) on all (right, wrong) pairs; report results."
    import time
    start = time.clock()
    good, unknown = 0, 0
    n = len(tests)
    for right, wrong in tests:
        w = correction(wrong)
        good += (w == right)
        if w != right:
            unknown += (right not in WORDS)
            if verbose:
                print('correction({}) => {} ({}); expected {} ({})'
                      .format(wrong, w, WORDS[w], right, WORDS[right]))
    dt = time.clock() - start
    print('{:.0%} of {} correct ({:.0%} unknown) at {:.0f} words per second '
          .format(good / n, n, unknown / n, n / dt))

spelltest(Testset(open('spell-testset1.txt')))  # Development set
spelltest(Testset(open('spell-testset2.txt')))  # Final test set
```

Measured output of the study:

```text
unit_tests pass
75% of 270 correct at 41 words per second     (development set)
68% of 400 correct at 35 words per second     (final test set)
```

Verbose-mode failures, the raw material for the next iteration (the counts
in parentheses are corpus frequencies - each line explains mechanically why
the wrong candidate won):

```text
correction('adres')  => 'acres' (37);   expected 'address' (77)
correction('rember') => 'member' (51);  expected 'remember' (162)
correction('thear')  => 'their' (3956); expected 'there' (4973)
correction('reciet') => 'recite' (5);   expected 'receipt' (14)
```

<!-- Reading the failure table by mechanism, not by instance: 'adres'->
'acres' loses because the trivial error model ranks ALL distance-1 edits
above ALL distance-2 edits, and doubling a letter (a common typo) is
distance 2 while d->c (a rare one) is distance 1 - so the fix is a weighted
error model, not more dictionary. 'thear' is undecidable at the single-word
level no matter the model - 'their' and 'there' are both one edit away and
both common - so that failure class argues for context (surrounding words),
a different axis entirely. Each future-work item is the shadow of an
observed failure class with a count attached; none is speculative. -->

## Distill - "Why Momentum Really Works" (Goh, 2017, CC-BY 4.0)

<!-- Annotation: The exposition register at its best. Three moves to copy.
FIRST: open on the reader's existing mental model and grade it - not "here is
our framework" but "here's the story you already believe, here's exactly what
it can't explain." SECOND: justify the model choice as minimal-sufficient.
THIRD: after proving the strongest possible claim (a matching lower bound),
immediately bound the claim's own scope. Attribution: Goh, "Why Momentum
Really Works", Distill, 2017, doi:10.23915/distill.00006, CC-BY 4.0. -->

The opening - the misconception, met head-on:

```text
Here's a popular story about momentum: gradient descent is a man walking down
a hill. He follows the steepest path downwards; his progress is slow, but
steady. Momentum is a heavy ball rolling down the same hill. The added
inertia acts both as a smoother and an accelerator, dampening oscillations
and causing us to barrel through narrow valleys, small humps and local
minima.

This standard story isn't wrong, but it fails to explain many important
behaviors of momentum. In fact, momentum can be understood far more precisely
if we study it on the right model.

One nice model is the convex quadratic. This model is rich enough to
reproduce momentum's local dynamics in real problems, and yet simple enough
to be understood in closed form.
```

The payoff, quantified rather than adjectival:

```text
The critical value of beta = (1 - sqrt(alpha*lambda_i))^2 gives us a
convergence rate (in eigenspace i) of 1 - sqrt(alpha*lambda_i). A square
root improvement over gradient descent, 1 - alpha*lambda_i!

[With both parameters optimal:]
    Convergence rate, Momentum:          (sqrt(k) - 1) / (sqrt(k) + 1)
    Convergence rate, Gradient Descent:  (k - 1) / (k + 1)

With barely a modicum of extra effort, we have essentially square rooted
the condition number!
```

The limits section - and the limit on the limit:

```text
Unfortunately, while improvements to the momentum algorithm do exist, they
all run into a certain, critical, almost inescapable lower bound. ...
[On Nesterov's worst-case function] the convergence rate that momentum
promises matches the best any linear first order algorithm can do. And we
arrive at the disappointing conclusion that on this problem, we cannot do
better.

Like many such lower bounds, this result must not be taken literally, but
spiritually. ... This lower bound does not preclude the possibility, for
example, of reformulating the problem to change the condition number itself!
There is still much room for speedups, if you understand the right places
to look.
```

<!-- Note the epistemics extend to the ARTICLE'S OWN REVIEW: Goh's published
acknowledgments credit an anonymous reviewer "for pointing out two
non-trivial errors in the original manuscript," with the review thread kept
public. Research writing that keeps its correction history visible is
modelling the same move as publishing the numbers that missed the target. -->

## Faraday - "The Chemical History of a Candle" (1861, public domain)

<!-- Annotation: Phenomenon-first, a century and a half before that was
advice. The reader watches a mundane object; Faraday extracts a mechanism
chain from what is literally visible: cup forms -> because the edge is cooler
-> because the rising current cools the sides -> because the flame heats the
air. Every link is checkable by the audience in the room. Note also "we
philosophers - I hope that I may class you and myself together in this case"
- the reader is enrolled as an investigator, not lectured at. -->

```text
You see that a beautiful cup is formed. As the air comes to the candle it
moves upwards by the force of current which the heat of the candle produces,
and it so cools all the sides of the wax, as to keep the edge much cooler
than the part within; the part within melts by the flame that runs down the
wick as far as it can go before it is extinguished, but the part on the
outside does not melt.
```

The lecture's stated program is a complete research agenda in one sentence:
to ascertain what happens in any particular part of the flame - why it
happens, what it does in happening, and where, after all, the whole candle
goes to. Observation, mechanism, dynamics, conservation - in that order.

---

## Negative pair (bad -> good) - the same study, two writeups

<!-- Original, written for this corpus, and FICTIONAL - continues the
multipart-upload universe of issues.md and review.md: after fixing the
signing bug, the team asks whether a larger chunk size is worth adopting.
The bad version is the AI-research-voice failure mode in concentrate: hedge
fog, no baseline, no numbers, criteria invented after the results, held-out
data silently tuned on, conclusions beyond the evidence. The good version is
the M&M/Norvig discipline applied to an ordinary engineering question. -->

**BAD** - hedge fog, no numbers, post-hoc criteria, grandiose conclusion:

```text
Title: Optimizing Upload Performance with Adaptive Chunking

It is well known that chunk size can significantly impact upload performance.
We explored various chunk sizes and it could be argued that the results may
potentially suggest that larger chunks provide substantially better
throughput in many cases. Our adaptive approach performed very well in our
testing and could revolutionize how uploads are handled. Memory usage was
acceptable. Future work includes exploring even more configurations.
```

**GOOD** - question and criteria first, method stated, exact results
including the failure, decision bounded by its limits:

```text
Question: does raising the multipart chunk size from 5 MB (current default)
to 16/32/64 MB improve wall-clock upload throughput, and at what memory cost?

Success criteria (declared before running anything): adopt the smallest chunk
size giving >= 10% median throughput gain on the 100 MB fixture, provided
peak RSS stays under 256 MB. Baseline: current 5 MB default.

Method: 30 uploads per configuration per network, interleaved round-robin to
absorb time-of-day drift. Two networks: office (1 Gbps fiber, ~0% loss) used
as the development environment, and one residential link (50 Mbps cable,
bursty loss) held out as the final test - untouched until the configuration
was chosen. foo-uploader 2.4.2, Node 20.11, eu-central-1.

Results, median MB/s over n=30 (IQR in brackets), 100 MB fixture:

               office 1 Gbps          home 50 Mbps (held out)
   5 MB        41.2  [39.8-42.6]      5.6  [5.3-5.8]     baseline
  16 MB        50.7  [49.1-52.0]      5.7  [5.4-5.9]     +23% / +2%
  32 MB        53.4  [50.2-55.1]      5.5  [4.9-5.8]     +30% / -2%
  64 MB        54.1  [48.7-56.3]      4.9  [3.1-5.6]     +31% / -12%

Peak RSS scaled linearly with chunk size as expected: 88 MB at 16 MB chunks,
301 MB at 64 MB (fails the memory criterion on its own).

Failure worth reporting: at 64 MB on the home link, 3/30 runs degenerated
into retry storms - a single lost segment forces retransmission of the whole
64 MB part, and the wide IQR (3.1-5.6) is that mechanism, not noise.

Decision: adopt 16 MB. It clears the 10% bar on the office network (+23%),
is throughput-neutral on the lossy link (+2%, within the IQR), and stays
well under the memory budget.

Limits (measured vs conjectured): measured on exactly two networks and one
OS; the loss-sensitivity mechanism predicts the 64 MB penalty GROWS with
loss rate, but we measured one lossy link, so that is a labelled conjecture,
not a result. RSS was sampled, not peak-committed memory. The 100 MB fixture
sits near the small end of production files (p50 = 340 MB); gains at p50 are
extrapolated, not measured.
```

<!-- Why the good one works: the criteria exist before the data, so the
decision is checkable rather than rationalized; the baseline and every
comparison carry exact magnitudes with spread; the held-out network plays
the role of the final test set and is labelled as such; the one ugly result
(retry storms) is reported WITH its mechanism, turning an anomaly into a
finding; and the limitations section separates what was measured from what
is being conjectured - M&M's certainly/probably distinction, and the same
fence issues.md puts around an unconfirmed hypothesis. The bad version
contains no sentence that could be false, which is exactly why it contains
no knowledge. -->