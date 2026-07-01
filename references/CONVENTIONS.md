# references/CONVENTIONS.md - paper-database spec

The only consumer of this database is an AI coding agent that **implements the papers' algorithms in C++/Python and cross-checks its code against them.**

Every rule below serves that and nothing else. Bibliographic completeness is a non-goal; fidelity to *implementable content* is the goal.

Read this file before ingesting a paper (the `paper-ingestor` subagent does) and before implementing from one (the consumer should).

---

## 1. Output layout (one folder per paper, split by PURPOSE not by section)

```text
references/
  index.yaml        # grep-first routing registry (all papers)
  papers/<key>/
    paper.md        # Tier 2: cleaned full body, ONE file, derivations intact
    equations.md    # Tier 1: implementable equations + co-located symbols
    notation.md     # Tier 1: symbol / meaning / units / conventions
    checks.md       # Tier 1: benchmark numbers -> "given X expect Y" fixtures
    meta.yaml       # this paper's index entry (also merged into index.yaml)
```

`<key>` = `firstauthor_year_topic`, lowercase, underscores (e.g. `sohldickstein_2014_lahmc`). Stable; used everywhere as a citation handle.

Do **not** split `paper.md` into per-section files - that severs derivations from their symbol scope. For very long book chapters, add an anchored `TOC` at the top of `paper.md` instead of fragmenting.

---

## 2. Keep / strip / extract taxonomy

| Element | Verdict | Destination |
|---|---|---|
| Equations defining quantities/algorithms | **extract** | `equations.md` (Tier 1) |
| Numeric results, result tables, test functions, hyperparameters, known constants | **extract** | `checks.md` - highest cross-check value |
| Symbols, units, sign/unit conventions (reduced units, `β=1/kT`, `mass=1`) | **extract** | `notation.md` |
| Method prose: assumptions, regime of validity, why | keep | `paper.md` |
| Derivations & proofs | keep, **label `Derivation`** | `paper.md` (NOT the equation sheet) |
| Figure captions carrying semantics ("Fig 1: F flips momentum") | reclassify to a note | `paper.md` / `notation.md` |
| Pure figure pointers ("see Fig 3") | strip | - |
| Title / authors / venue / year / DOI / ISBN | extract | `meta.yaml` |
| References / bibliography | strip body, **harvest keys** first | `depends_on` in `meta.yaml` |
| Intro literature-review citation soup | strip prose, harvest keys | `depends_on` |
| Affiliations, emails, dates, correspondence, copyright, funding, acknowledgments | strip | raw only |
| Page numbers, running heads, `<span id=...>`, `![](_page_...)` image refs | strip | raw only |

---

## 3. The two paper species (classify first)

- **methods / theory paper** (rich implementable math): full Tier-1 extraction.
- **software / survey / overview paper** with **no implementable equations** (e.g. a library-description paper whose value is pointing at *other* papers):
  - set `has_implementable_math: false`, **do not fabricate an equation sheet**,
  - emit only `paper.md` + `meta.yaml` (+ `checks.md` if it states any numbers).
  - Its value is routing: fill `depends_on` and `provides` well.

If you find yourself inventing equations to fill `equations.md`, stop - the paper is probably routing-only.

---

## 4. Grep / LLM-friendliness rules (apply while cleaning)

- Normalize display math to `$$ ... $$`, inline to `$ ... $`. Convert `\begin{equation}` and `\[ ... \]`; leave `\( \)` alone (this converter uses them for escaped parens in citations, not math).
- Above every display equation write an anchor comment: `<!-- eq:LABEL -->` (use the paper's own number if present: `eq:25`). Grep the anchor, not the raw TeX.
- **Co-locate symbol definitions with each equation.** An equation without its symbols is uncodeable.
- Real `##` headers per section. Wrap prose to paragraphs, not 5000-char lines.

### OCR corruption (important for these files)

- Display `$$...$$` blocks usually survive conversion; **inline math is often corrupted** - `T (x 0 |x)` means `T(x'|x)`, `0` is frequently a lost prime, integrals/sums drop out.
- Repair using surrounding prose, and mark every uncertain fix with `<!-- CHECK: ... -->`.
- "Renders valid" != "semantically correct": a clean-parsing equation can still be wrong, so sanity-check the key ones against the prose, not just the brace-checker.

---

## 5. Validation

- Run `python references/tools/check_latex.py <file.md>` after writing `equations.md`. It reports unbalanced braces / `\left`\`\right` / `\begin`\`\end`. Fix or explicitly`<!-- CHECK -->`-flag every failure before finishing.

---

## 6. Dedup (check BEFORE doing expensive extraction)

A paper already exists if any of, in order:

1. same normalized DOI,
2. fuzzy title match ≥ 0.90 **and** (same year **or** title match ≥ 0.94)
   - harvested year is unreliable, so a very strong title match wins alone.
If it exists and `--force` was not requested, skip and report.

---

## 7. Book chapters

- No abstract / no DOI / "Chapter N" structure -> `doc_type: book_chapter`.
- Use `isbn`, `book_title`, `chapter` instead of `doi`/`venue`. Chapters are derivation-heavy (bigger `Derivation` blocks) and long (anchored `TOC`, don't fragment).
- Exercises -> keep in `paper.md` as `Exercises`; treat as checks only if answers are given.

---

## 8. index.yaml entry schema (capability-first, NOT bibliography-first)

```yaml
- key: sohldickstein_2014_lahmc
  doc_type: methods            # article|methods|survey|software|book_chapter|thesis
  title: "Hamiltonian Monte Carlo Without Detailed Balance"
  authors: [Sohl-Dickstein, Mudigonda, DeWeese]
  year: "2014"
  venue: ICML
  doi: null
  summary: "HMC variant that avoids rejection by extending trajectories; drop-in for samplers."
  has_implementable_math: true
  methods: [leapfrog integrator, momentum flip, Look-Ahead HMC transition]
  provides: [LAHMC sampler, transition-probability rule, reduced momentum flips]
  key_equations: [{label: "23", desc: "LAHMC transition operator"},
                  {label: "25", desc: "greedy leapfrog transition probabilities"}]
  keywords: [MCMC, HMC, detailed balance, sampling]
  domains: [molecular-sim, statistics]
  depends_on: [duane_1987, neal_2010, hairer_2003]   # routing edges
  path: references/papers/sohldickstein_2014_lahmc
  content_hash: "..."
```

`methods` / `provides` / `key_equations` are the retrieval surface - the consumer searches by *capability* ("I need a leapfrog integrator"), not by DOI.
