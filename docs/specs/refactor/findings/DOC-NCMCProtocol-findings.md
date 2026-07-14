# DOC-NCMCProtocol findings

## Refined hypothesis (imprecise, not contradicted)
- Ticket section 1 asked to document "lambda(0), lambda(N), and monotonicity as
  contract." The schedule is NOT globally monotonic: it is a PALINDROMIC tent
  (`1 -> 0 -> 1`), monotone decreasing to the center then monotone increasing.
  Verified against `tests/TestNCMCWork.cpp`:
  - endpoints pinned: `protocolLambda(0,n,.)==protocolLambda(n-1,n,.)==1`
    (`TestNCMCWork.cpp:81,93,94`);
  - palindrome: `protocolLambda(s)==protocolLambda(n-1-s)` (`:109`);
  - piecewise-monotone down (`:150`) then up (`:153`), center value 0 (`:136`);
  - range `[0,1]` (`:122,183`); degenerate `protocolLambda(0,1,.)==1` (`:179`).
  Documented the palindrome + pinned endpoints as the contract, and stated the
  monotonicity precisely (piecewise, symmetric) rather than "monotonic".

## Endpoint convention (verified)
- `lambda == 1` fully coupled, `lambda == 0` fully decoupled; both NCMC endpoints
  are fully coupled. Caller: `src/world/sampler/NcmcMove.cpp:72-77` (World
  delegates) and `:381` (per-substep use). The alchemy force SHALL share this
  convention; documented as a `@file` invariant.

## Notes
- No `@note Assumed:` used.
