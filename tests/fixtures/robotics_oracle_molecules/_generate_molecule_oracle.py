"""Clone-side Scope-B molecule oracle generator (docs/specs/robotics-oracle-
differential.md §4.2 Scope-B build path, resolved follow-up 3).

For each of the three §4.4 molecular topology classes:
  * rigid   -- examples/10ala.*, all bonds Rigid (Weld root; the port's
               companion test is RoboticsOracleMolecule.RigidWeldRoot).
  * regular -- examples/10ala.*, default-flexible (every non-terminal,
               non-ring bond -> Torsion), matching the PORT's own
               ``Context::buildFlexibilities(nullopt, Torsion, false)``
               selection via ``context.selectBonds("all", excludeTerminal=True)``
               (both pick every bond whose two endpoints are non-terminal;
               ring-closing bonds are forced back to Rigid by the clone's own
               ``Context::addWorld``, see src/Context.cpp, matching the port).
  * cyclic  -- examples/1APQ.*, same flexible selection as "regular".

this builds the molecule through the clone's OWN Python pipeline
(``Context.add_robotic_world`` / ``selectBonds`` / ``create_torsional_bonds`` --
exactly what ``add_robotic_world``'s docstring and ``test_installation.py``
exercise), then calls the additive
``World.dump_robotics_oracle_molecule`` (src/PyBind11.cpp ->
Robosample/src/RoboticsOracleMoleculeDump.cpp) for the "rest" state and (for
the regular/cyclic classes, which have >=1 Torsion dof) one seeded "random"
state, and writes ``<case>.moldyn.npz`` + ``<case>.moldyn.manifest.json``
under this directory -- alongside the port-side ``<mol>.systopo.npz`` fixtures
``_generate_port_topology.py`` already produces, which the port comparator
(tests/TestRoboticsOracleMolecule.cpp) loads via ``World::buildModel`` to run
the SAME differential.

The rigid (Weld-root) class has zero internal DOF, so a "random" state is not
meaningful there (see docs/specs/robotics-oracle-differential.md: the
RigidFreeRoot numeric differential is deliberately out of scope for THIS
generator too -- the clone's simplified ``Context.add_robotic_world`` exposes
no root-mobility knob, so a Free-root molecule build would require the
separate ``addDockingWorld`` path; the structural Free-root check already
lives in the port-only ``RoboticsOracleMolecule.RigidFreeRoot`` test).

Manual dev tool (no ``test_`` prefix -> never pytest-collected). Run by hand:

    python3 tests/fixtures/robotics_oracle_molecules/_generate_molecule_oracle.py

Requires the CLONE's ``robo_bindings`` extension to be built AND installed
into ``Robosample/python/robosample/`` (``cmake --build --preset
cuda-release`` then ``cmake --install build/cuda-release``, run from
``Robosample/``, per docs/specs/robotics-oracle-differential.md §11 step 1).
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]  # .../Robosample (this, disasm repo root)
CLONE_ROOT = REPO_ROOT / "Robosample"  # .../Robosample/Robosample (the "refactor" clone)
FIXTURE_DIR = pathlib.Path(__file__).resolve().parent

# The CLONE's python package, NOT the port's (REPO_ROOT / "python"), which
# _generate_port_topology.py uses -- both are importable as `robosample`, so
# this script must run in its own fresh interpreter (never import both in one
# process; each generator is a separate `python3 <script>.py` invocation).
sys.path.insert(0, str(CLONE_ROOT / "python"))
import robosample as rb  # noqa: E402

# §8.3 convention: no live RNG at gate time from the PORT's perspective -- the
# clone (which owns this oracle) draws the randoms ONCE, deterministically
# (fixed seed), and bakes both the inputs (q, u) and the outputs into the
# fixture; the port replays the baked q/u verbatim (see
# RoboticsOracleMoleculeDump.hpp's "random" state contract).
RNG_SEED = 20260702

MOLECULES = {
    "10ala": (REPO_ROOT / "examples" / "10ala.prmtop", REPO_ROOT / "examples" / "10ala.rst7"),
    "1APQ": (REPO_ROOT / "examples" / "1APQ.prmtop", REPO_ROOT / "examples" / "1APQ.rst7"),
}


def _load_port_ring_closing_prmtop_pairs(mol_name: str) -> list[tuple[int, int]]:
    """Read the PORT's own ring-closing bond set (docs/specs/robotics-oracle-
    differential.md Scope B §6 shared-tree fix) out of the ``<mol>.systopo.npz``
    fixture ``_generate_port_topology.py`` already wrote, keyed by raw prmtop
    atom index (the cross-engine identity currency, §5) -- NOT the port's own
    internal (BFS-order) atom index that ``bonds_i``/``bonds_j`` are stored in.

    Plain ``numpy.load`` only -- no ``robosample`` import of either flavor is
    needed to read a ``.npz``, so this stays safe to call from this script
    (which has the CLONE's ``robosample`` on ``sys.path``, never the port's).
    """
    npz = np.load(FIXTURE_DIR / f"{mol_name}.systopo.npz")
    bonds_i = npz["bonds_i"]
    bonds_j = npz["bonds_j"]
    ring_closing = npz["bonds_ring_closing"]
    prmtop_index = npz["atoms_prmtop_index"]
    pairs = [
        (int(prmtop_index[bonds_i[k]]), int(prmtop_index[bonds_j[k]]))
        for k in range(len(ring_closing))
        if ring_closing[k]
    ]
    return pairs


def _new_context(
    case_name: str,
    prmtop: pathlib.Path,
    inpcrd: pathlib.Path,
    ring_closing_bond_prmtop_pairs: list[tuple[int, int]] | None = None,
) -> rb.Context:
    return rb.Context(
        name=case_name,
        seed=0,
        prmtop=str(prmtop),
        inpcrd=str(inpcrd),
        write_freq=1,
        testing=False,
        ring_closing_bond_prmtop_pairs=ring_closing_bond_prmtop_pairs,
    )


def _assert_ring_closing_bonds_match_port(context: "rb.Context", port_pairs: list[tuple[int, int]]) -> None:
    """§1 of the shared-tree fix: a REAL set-equality check (not a count),
    keyed by prmtopIndex bond endpoints, that this (clone) engine's ACTUAL
    ring-closing bond set -- as accepted by ``MoleculePrototype``/Molmodel,
    read back from ``context.system_topology.bonds`` -- matches the port's.
    This is the load-bearing assertion: it is what proves the override in
    ``ring_closing_bond_prmtop_pairs`` above was not silently dropped or
    partially rejected (e.g. by Molmodel's own bond-center bookkeeping),
    not merely that we asked for the right set.
    """
    port_set = {frozenset(p) for p in port_pairs}
    clone_set = {
        frozenset(b.prmtop_indices) for b in context.system_topology.bonds if b.ring_closing
    }
    if clone_set != port_set:
        missing_in_clone = port_set - clone_set
        extra_in_clone = clone_set - port_set
        raise AssertionError(
            "clone ring-closing bond set does NOT match the port's (shared-tree "
            f"fix failed): missing_in_clone={sorted(map(sorted, missing_in_clone))} "
            f"extra_in_clone={sorted(map(sorted, extra_in_clone))}"
        )


def _finalize_world(context: rb.Context) -> "rb.World":
    """add_sampler + initialize is required by Context.initialize (asserts
    exactly one sampler per world) even though this dump never runs any
    sampling -- mdSteps=0 is a harmless placeholder."""
    context.getWorld(0)  # no-op; keeps the intent explicit that world 0 is what we return below
    return context.getWorld(0)


def dump_case(
    case_name: str,
    molecule_class: str,
    context: rb.Context,
    want_random: bool,
    ring_closing_bond_prmtop_pairs: list[tuple[int, int]] | None = None,
) -> None:
    world = context.getWorld(0)

    states = ["rest"]
    if want_random:
        states.append("random")

    dumps = {label: world.dump_robotics_oracle_molecule(label, RNG_SEED) for label in states}

    # Body layout (atom-set, nq, nu) is constant across states within a case
    # -- take it from the first state and sanity-check it does not drift.
    ref_bodies = dumps[states[0]]["bodies"]
    for label in states[1:]:
        bodies = dumps[label]["bodies"]
        assert len(bodies) == len(ref_bodies), f"{case_name}: body count changed between states"
        for b_ref, b in zip(ref_bodies, bodies):
            assert b_ref["nq"] == b["nq"] and b_ref["nu"] == b["nu"], (
                f"{case_name}: per-body nq/nu changed between states -- not the same molecule build"
            )
            assert b_ref["atom_prmtop_indices"] == b["atom_prmtop_indices"], (
                f"{case_name}: per-body atom set changed between states"
            )

    npz_arrays: dict[str, np.ndarray] = {}
    for s, label in enumerate(states):
        d = dumps[label]
        npz_arrays[f"state{s}_logDetM"] = np.asarray([d["log_det_m"]], dtype=np.float64)
        npz_arrays[f"state{s}_KE"] = np.asarray([d["kinetic_energy"]], dtype=np.float64)
        npz_arrays[f"state{s}_normUdot"] = np.asarray([d["norm_udot"]], dtype=np.float64)
        npz_arrays[f"state{s}_minEigD"] = np.asarray([d["min_eig_d"]], dtype=np.float64)
        for b, bd in enumerate(d["bodies"]):
            prefix = f"state{s}_body{b}_"
            npz_arrays[prefix + "q"] = np.asarray(bd["q"], dtype=np.float64)
            npz_arrays[prefix + "u"] = np.asarray(bd["u"], dtype=np.float64)
            npz_arrays[prefix + "udot"] = np.asarray(bd["udot"], dtype=np.float64)
            # Frame-invariant per-atom Ground position anchor (docs/specs/
            # robotics-oracle-differential.md refined §6): flattened xyz,
            # PARALLEL to atom_prmtop_indices (constant across states, stored
            # once in the manifest's "bodies" list below).
            npz_arrays[prefix + "atomPosG"] = np.asarray(bd["atom_pos_g"], dtype=np.float64)
            npz_arrays[prefix + "X_GB_R"] = np.asarray(bd["X_GB_R"], dtype=np.float64)
            npz_arrays[prefix + "X_GB_p"] = np.asarray(bd["X_GB_p"], dtype=np.float64)
            npz_arrays[prefix + "V_GB_ang"] = np.asarray(bd["V_GB_ang"], dtype=np.float64)
            npz_arrays[prefix + "V_GB_lin"] = np.asarray(bd["V_GB_lin"], dtype=np.float64)
            npz_arrays[prefix + "A_GB_ang"] = np.asarray(bd["A_GB_ang"], dtype=np.float64)
            npz_arrays[prefix + "A_GB_lin"] = np.asarray(bd["A_GB_lin"], dtype=np.float64)
            # Zero-force mobilizer reaction (docs/specs/robotics-oracle-
            # reactions.md §4): the quantity World::calcSpatialForces stores,
            # reconstructed at zero applied force. Bo = convention-free
            # cross-engine anchor; Mo.angular = convention-gated. schema v2.
            npz_arrays[prefix + "reactionBo_ang"] = np.asarray(bd["reaction_bo_ang"], dtype=np.float64)
            npz_arrays[prefix + "reactionBo_lin"] = np.asarray(bd["reaction_bo_lin"], dtype=np.float64)
            npz_arrays[prefix + "reactionMo_ang"] = np.asarray(bd["reaction_mo_ang"], dtype=np.float64)
            npz_arrays[prefix + "reactionMo_lin"] = np.asarray(bd["reaction_mo_lin"], dtype=np.float64)
            # Per-body min-eig(D_b) (docs/specs/singular-dof-fixman.md): lets
            # the port-side test guard the axis-direction-dependent GAUGE
            # comparisons (V_GB.angular / A_GB.angular) per-body, rather than
            # tree-wide via the aggregate state{s}_minEigD above.
            npz_arrays[prefix + "minEigD"] = np.asarray([bd["min_eig_d"]], dtype=np.float64)

    npz_path = FIXTURE_DIR / f"{case_name}.moldyn.npz"
    np.savez(npz_path, **npz_arrays)

    # This engine's ACTUAL ring-closing (cotree) bond set, read back from
    # ``system_topology.bonds`` (authoritative -- reflects what Molmodel
    # actually accepted, not merely what was requested via
    # ``ring_closing_bond_prmtop_pairs``). Baked into the manifest so the
    # port-side C++ test (tests/TestRoboticsOracleMolecule.cpp) can assert
    # SET equality (keyed by prmtopIndex bond endpoints) against its own
    # ``SystemTopology.bonds_ring_closing`` on every ``ctest`` run -- the
    # real, durable form of the shared-tree fix's §1 check (docs/specs/
    # robotics-oracle-differential.md Scope B §6), not just a one-off
    # generation-time sanity check.
    clone_ring_closing_pairs = sorted(
        sorted(b.prmtop_indices) for b in context.system_topology.bonds if b.ring_closing
    )
    if ring_closing_bond_prmtop_pairs is not None:
        _assert_ring_closing_bonds_match_port(context, ring_closing_bond_prmtop_pairs)

    manifest = {
        "schema_version": 2,
        "case_name": case_name,
        "molecule_class": molecule_class,  # "rigid" | "regular" | "cyclic"
        "rng_seed": RNG_SEED,
        "states": states,
        "num_bodies": len(ref_bodies),
        "ring_closing_bond_prmtop_pairs": clone_ring_closing_pairs,
        "bodies": [
            {
                "index": b,  # 0-based; MobilizedBodyIndex == index + 1
                "nq": bd["nq"],
                "nu": bd["nu"],
                "atom_prmtop_indices": bd["atom_prmtop_indices"],
            }
            for b, bd in enumerate(ref_bodies)
        ],
        "arrays": sorted(npz_arrays.keys()),
    }
    manifest_path = FIXTURE_DIR / f"{case_name}.moldyn.manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2))

    print(
        f"wrote {npz_path} + {manifest_path} "
        f"(class={molecule_class}, numBodies={len(ref_bodies)}, states={states})"
    )


def main() -> None:
    # ---- Rigid: 10ala, all-Rigid (Weld root; 0 DOF -> rest state only) -----
    prmtop, inpcrd = MOLECULES["10ala"]
    ctx = _new_context("10ala_rigid", prmtop, inpcrd)
    ctx.add_robotic_world([[]]).add_sampler(
        timeStep=0.001, mdSteps=0, boostMDSteps=0,
        acceptRejectMode=rb.rb.AcceptRejectMode.MetropolisHastings,
    )
    ctx.initialize([300.0])
    dump_case("10ala_rigid", "rigid", ctx, want_random=False)

    # ---- Regular: 10ala, default-flexible (Torsion) ------------------------
    prmtop, inpcrd = MOLECULES["10ala"]
    ctx = _new_context("10ala_regular", prmtop, inpcrd)
    flex = [ctx.selectBonds("all", excludeTerminal=True)]
    ctx.add_robotic_world(flex).add_sampler(
        timeStep=0.001, mdSteps=0, boostMDSteps=0,
        acceptRejectMode=rb.rb.AcceptRejectMode.MetropolisHastings,
    )
    ctx.initialize([300.0])
    dump_case("10ala_regular", "regular", ctx, want_random=True)

    # ---- Cyclic: 1APQ, default-flexible (Torsion; ring bonds forced Rigid) -
    # Shared-tree fix (docs/specs/robotics-oracle-differential.md Scope B §6,
    # hostile-review verdict): 1APQ has 13 ring-closing bonds (prolines +
    # disulfide CYX cross-links), and this engine's OWN cycle-basis/maximum-
    # spanning-tree heuristic may pick a DIFFERENT (also topologically valid)
    # spanning tree than the port's amber_loader.py did -- making a per-body
    # q/u copy across engines invalid at q!=0 even though both trees are
    # individually correct. Drive this engine to break the SAME bonds the
    # port already broke (read back from the port's OWN dumped
    # SystemTopology, tests/fixtures/robotics_oracle_molecules/
    # _generate_port_topology.py's 1APQ.systopo.npz) so the two engines'
    # spanning trees are IDENTICAL, not merely both valid.
    prmtop, inpcrd = MOLECULES["1APQ"]
    port_ring_closing_pairs = _load_port_ring_closing_prmtop_pairs("1APQ")
    ctx = _new_context(
        "1APQ_cyclic", prmtop, inpcrd,
        ring_closing_bond_prmtop_pairs=port_ring_closing_pairs,
    )
    flex = [ctx.selectBonds("all", excludeTerminal=True)]
    ctx.add_robotic_world(flex).add_sampler(
        timeStep=0.001, mdSteps=0, boostMDSteps=0,
        acceptRejectMode=rb.rb.AcceptRejectMode.MetropolisHastings,
    )
    ctx.initialize([300.0])
    dump_case(
        "1APQ_cyclic", "cyclic", ctx, want_random=True,
        ring_closing_bond_prmtop_pairs=port_ring_closing_pairs,
    )


if __name__ == "__main__":
    main()
