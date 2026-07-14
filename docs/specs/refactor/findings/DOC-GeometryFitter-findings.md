# DOC-GeometryFitter findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed
in GeometryFitter.cpp.

## Edit
The full INV-3 reconstruction contract lived at the definition
(GeometryFitter.cpp:89-141) but the public declaration in World.hpp was bare.
Added a concise contract at the `setAtomsLocationsInGround` /
`getAtomsLocationsInGround` declarations (INV-3 reconstruct-from-coordinates,
borrow/lifetime, Cartesian-branch behavior). The definition narrative (Gibbs
continuation rationale) was left intact; the header carries the crisp contract,
the definition the rationale -- complementary, not duplicated.

## Verified hypotheses
- INV-3 sole-currency + full reconstruction: verified. q reset to 0 (identity
  quaternion) does not discard geometry; frames are rebuilt from the incoming
  Cartesian, so the reset only re-references the block (GeometryFitter.cpp
  :100-110). Round-trip set->get is identity up to the tree fit.
- Loop-closure targets track the carried-over distance, not force-field r0
  (:130-140) -- verified; part of the continued state.

## Assumed notes
None.
