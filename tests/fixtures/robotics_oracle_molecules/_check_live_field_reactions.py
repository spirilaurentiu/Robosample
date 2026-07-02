"""Clone-side Leg B check: live-force-field mobilizer reaction validation
(docs/specs/robotics-oracle-reactions.md §5).

Validates exactly the quantity ``World::calcSpatialForces`` stores -- the
mobilizer reaction ``matter->calcMobilizerReactionForces`` under the LIVE
DuMM/OpenMM force field on a REAL molecule -- against Simbody's INDEPENDENT
``calcMobilizerReactionForcesUsingFreebodyMethod`` (a from-scratch Newton-Euler
freebody sweep, ~3x slower, structurally independent of the fast PPlus path).
Unlike the Scope-B differential (Leg A, zero applied force, §4.5), this
exercises the per-atom->per-body force reduction under a live field.

This is the ONE reaction check that runs under a live force field, so it is
the direct guardian of ``calcSpatialForces``'s stored ``force_G``/``torque_G``.

Manual dev tool (no ``test_`` prefix -> never pytest-collected), CLONE-side:
it needs the clone's ``robo_bindings`` extension built + installed (the disasm
``nox -s tests`` gate does NOT build clone bindings, docs/specs/robotics-oracle-
reactions.md §5.2 NOTE), so it lives with ``_generate_molecule_oracle.py`` and
is run by hand after a clone build:

    python3 tests/fixtures/robotics_oracle_molecules/_check_live_field_reactions.py

Exits nonzero (raises AssertionError) if any molecule's two-method residual
exceeds tolerance, so it can be wired into a clone-side gate later.
"""

from __future__ import annotations

import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]  # .../Robosample (disasm repo root)
CLONE_ROOT = REPO_ROOT / "Robosample"  # .../Robosample/Robosample (the "refactor" clone)

# The CLONE's python package (see _generate_molecule_oracle.py's note): run in
# a fresh interpreter, never import both flavors of `robosample` in one process.
sys.path.insert(0, str(CLONE_ROOT / "python"))
import robosample as rb  # noqa: E402

RNG_SEED = 20260702  # unused here (no random state), kept for parity with the generator

MOLECULES = {
    "10ala": (REPO_ROOT / "examples" / "10ala.prmtop", REPO_ROOT / "examples" / "10ala.rst7"),
    "1APQ": (REPO_ROOT / "examples" / "1APQ.prmtop", REPO_ROOT / "examples" / "1APQ.rst7"),
}

# The two Simbody methods are documented to "agree to within numerical
# precision" (SimbodyMatterSubsystem.h:2940). On a real molecule the freebody
# tip-to-base sweep accumulates O(depth) round-off, so use a modest absolute
# floor on the relative residual -- still ~7 orders below an O(1) transcription
# bug in the reaction path.
REL_TOL = 1e-9


def _new_context(name: str, prmtop: pathlib.Path, inpcrd: pathlib.Path) -> rb.Context:
    return rb.Context(
        name=name, seed=0, prmtop=str(prmtop), inpcrd=str(inpcrd),
        write_freq=1, testing=False, ring_closing_bond_prmtop_pairs=None,
    )


def check_case(name: str, prmtop: pathlib.Path, inpcrd: pathlib.Path) -> None:
    ctx = _new_context(name, prmtop, inpcrd)
    flex = [ctx.selectBonds("all", excludeTerminal=True)]
    ctx.add_robotic_world(flex).add_sampler(
        timeStep=0.001, mdSteps=0, boostMDSteps=0,
        acceptRejectMode=rb.rb.AcceptRejectMode.MetropolisHastings,
    )
    ctx.initialize([300.0])
    world = ctx.getWorld(0)

    res = world.check_live_field_reactions()
    print(
        f"[{name}] numBodies={res['num_bodies']} "
        f"maxAbsDiff={res['max_two_method_abs_diff']:.3e} "
        f"maxRelDiff={res['max_two_method_rel_diff']:.3e} "
        f"maxReactionNorm={res['max_reaction_norm']:.3e} "
        f"worstBody={res['worst_body']}"
    )
    assert res["max_two_method_rel_diff"] <= REL_TOL, (
        f"[{name}] live-field reaction FAILED: the fast calcMobilizerReactionForces "
        f"(what World::calcSpatialForces stores) disagrees with the independent "
        f"freebody method by rel={res['max_two_method_rel_diff']:.3e} > {REL_TOL:.0e} "
        f"at body {res['worst_body']} -- a real reaction-path bug under the live force field"
    )
    # A live force field on a flexible tree must transmit a NONZERO reaction
    # somewhere (else the check is vacuous, Rule 8).
    assert res["max_reaction_norm"] > 1.0, (
        f"[{name}] max reaction norm {res['max_reaction_norm']:.3e} is implausibly small -- "
        f"the live force field is not reaching the reactions (vacuous check)"
    )


def main() -> None:
    check_case("10ala_live", *MOLECULES["10ala"])
    check_case("1APQ_live", *MOLECULES["1APQ"])
    print("live-field reaction check PASSED")


if __name__ == "__main__":
    main()
