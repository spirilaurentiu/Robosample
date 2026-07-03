# NCMC explicit-solvent acceptance — Region and protocol policy

## 1. Scope
Region-selection and electrostatics-decoupling policy for a **flexible interesting molecule making an internal torsional basin-hop** (not a small rigid ligand). Answers mission Q5. Efficiency, not correctness (all options are exact at `λ=1`). Terms from `00-diagnosis-and-scaling.md` §2.

## 2. Region A: decouple the moving substructure, not the whole molecule
**POLICY (NORMATIVE for flexible receptors).** Region A SHALL be the substructure whose intermolecular contacts actually clash during the target torsional move (the loop/sidechain undergoing the hop plus its immediate intermolecular contacts), not the entire receptor.

Rationale. `⟨W_diss⟩` scales with the solvation-shell area of A that must re-solvate during the switch (`10-acceptance-construction.md` §4). Decoupling the whole receptor re-solvates the entire receptor surface every move — enormous, mostly wasted reorganization work, since most of the receptor is not moving. A smaller A means a smaller shell, fewer stages for a reversible switch, and lower acceptance cost.

Constraint. A need not be a whole molecule. **DECIDED: the alchemy API SHALL accept an arbitrary atom-INDEX SET, not only a contiguous `[begin,end)` block** — `configure_ncmc`/`add_ncmc_world` (`context.py:948–960`) and `createAlchemyDecouplingForces` (`OpenMMContext.cpp:725`) are extended so Region A is any index set, letting a non-contiguous moving substructure (loop/sidechain) be decoupled. The soft-core force already takes arbitrary sets via `addInteractionGroup(aSet, restSet)` (`:762`); the PME charge-offset loop (`:738`) and the intra-A hard-LJ set (`:759–768`) SHALL iterate the index set rather than `[begin,end)`. Bonded terms crossing A's boundary stay physical; only *intermolecular* nonbonded to B is scaled. NOTE: preserve a contiguous-range convenience overload for the existing callers.

## 3. Intra-region electrostatics: decouple, do not annihilate
**Current behavior.** `createAlchemyDecouplingForces` (`OpenMMContext.cpp:725`) scales A's PME charge via a parameter offset (`:742`) so charge(λ)=λ·q *against everything*, which **annihilates intra-A electrostatics for λ<1** (documented at `:729–734`). LJ is handled decoupling-correctly (soft-core A×rest at `:752`, hard intra-A LJ restored at `:768`), but electrostatics is annihilation, not decoupling.

**CLAIM P1.** For a flexible A this distorts the trough proposal and hurts acceptance (efficiency, not correctness). During the λ=0 uncaged stride A's *internal* electrostatics vanish, so the torsional landscape sampled at the trough is missing the intramolecular Coulomb (i,i+4 H-bonds, salt bridges, rotamer-shaping terms) that define the real internal basins. The stride proposes internal states favorable *without* intramolecular electrostatics; those are then rejected at recoupling. The larger and more polar A is, the worse the mismatch.

**POLICY (RECOMMENDED).** For a flexible A the intra-A electrostatics SHALL be preserved across λ (intermolecular-only decoupling): keep A's charges physical at all λ and scale only the A↔B Coulomb. Implementation options for the reviewer:
- add a `CustomNonbondedForce`/`CustomBondForce` intra-A Coulomb (with A's exclusions) that is λ-independent, mirroring the hard intra-A LJ pattern (`:768`), and restrict the PME charge-offset so it scales only A↔B; or
- a reaction-field/decoupling electrostatics scheme where the reciprocal-space A self-term is compensated.

Annihilation is acceptable only for a small rigid ligand whose internal landscape is trivial (no distortion to hurt). NOTE this is efficiency policy; both schemes are exact at `λ=1`.

## 4. Teleport is a small-ligand device, out of scope here
`ncmcApplyTroughTeleport` (`:1842`) rigidly repositions a Free-root region at the λ=0 trough — for a *small mobile molecule* with a Free root and a defined docking site. A flexible receptor doing an internal torsional hop has its root welded to Ground and no Free root to teleport; the basin-hop is driven by the torsional DOF integrating freely during the λ=0 hold. **POLICY:** teleport SHALL remain inactive for the flexible-region use (the existing guard `constraints_.empty() && Free-root && siteAtoms` already disables it); it is not part of this fix.

## 5. Touch list
- `src/OpenMMContext.cpp:725–800` (`createAlchemyDecouplingForces`), `:617` (`setAlchemicalLambda`) — intra-A electrostatics preservation. Conventions at risk: PME reciprocal-space correctness of the charge-offset route; the shared exclusion list (`addStandardExclusions`) must still make intra-A pairs non-double-counted.
- `python/robosample/context.py add_ncmc_world` — region selection to a moving-substructure block.
