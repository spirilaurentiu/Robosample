"""Tier 2 (docs/specs/fixman-idealized-chains-validation.md): proline
dipeptide realism -- vacuum AMBER ff on Ace-Pro-Pro-Nme (2 loop-closure
constraints; the spec's "richer valid variant", the files already shipped
under ``examples/fixman/``). Python full-stack path (spec Sec. 3.3).

T2.0 is the always-on structural guard (no slow-tier gate): assert the cyclic
Fixman path (``World::numLoopConstraints`` / ``Constraints::
calcConstraintLogDet``, exposed to Python as ``World.num_loop_constraints`` /
``World.current_constraint_log_det()`` -- a small, surgical PyBind11/World.hpp
addition; see the coder checkpoint) actually fires and is conformation-
dependent -- neither was reachable from Python before, so without this guard
the cyclic path could silently go untested. T2.1 (slow) is the full
statistical Hellinger comparison, gated by ``ROBOSAMPLE_SLOW_TESTS`` and run
at a REDUCED sample count for the same per-round-cost reason as
``test_fixman_idealized_chains.py`` (see that module's docstring).
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

pytest.importorskip("openmm")
robosample = pytest.importorskip("robosample")
pytest.importorskip("mdtraj")

from fixman_chain_lib import (  # noqa: E402
    build_cartesian_context,
    build_torsional_context,
    hellinger_distance,
    histogram_counts,
    run_and_load_dihedrals,
    slow_tests_enabled,
    thin_to_independent,
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PRMTOP = REPO_ROOT / "examples" / "fixman" / "ace_pro_pro_nme.prmtop"
RST7 = REPO_ROOT / "examples" / "fixman" / "ace_pro_pro_nme.rst7"


def _bfs_adjacency(context):
    st = context.system_topology
    adj: dict[int, list[int]] = {}
    for i, j in zip(list(st.bonds_i), list(st.bonds_j)):
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    return adj


def _dihedral_quads_prmtop_order(context, bonds):
    """For each rotatable central bond (i,j) (BFS/compound atom-index order --
    RobotModel.hpp's own index convention, NOT prmtop order), complete it to a
    dihedral quad (h,i,j,k) using any other bonded neighbor on each side, then
    remap every index through ``atoms_prmtop_index`` -- the DCD trajectory
    (and hence mdtraj's topology, loaded from the .prmtop) is written in
    PRMTOP atom order (``Context::writeOutputs``' ``atomsPrmtopIndex``
    scatter), which differs from the BFS/compound order for a branched,
    multi-residue molecule like a peptide (unlike the purely-linear bead
    chains in test_fixman_idealized_chains.py, where the two orders happen to
    coincide).
    """
    st = context.system_topology
    prmtop_idx = list(st.atoms_prmtop_index)
    adj = _bfs_adjacency(context)
    quads = []
    for i, j in bonds:
        h = next((n for n in adj.get(i, []) if n != j), None)
        k = next((n for n in adj.get(j, []) if n != i), None)
        if h is None or k is None:
            continue
        quads.append([prmtop_idx[h], prmtop_idx[i], prmtop_idx[j], prmtop_idx[k]])
    return quads


def _skip_if_missing():
    if not PRMTOP.exists() or not RST7.exists():
        pytest.skip(f"example not found: {PRMTOP}")


# ---------------------------------------------------------------------------
#  T2.0 -- cyclic-path structural guard (always-on). SHALL.
# ---------------------------------------------------------------------------
def test_ace_pro_pro_nme_loop_closure_structural_guard(tmp_path):
    _skip_if_missing()

    ctx, world, bonds = build_torsional_context(
        str(tmp_path / "pro_t20"),
        PRMTOP,
        RST7,
        seed=101,
        use_fixman=True,
        timestep=0.0007,
        md_steps=8,
    )
    assert len(bonds) >= 1
    ctx.initialize([300])

    # Ace-Pro-Pro-Nme is DI-proline: two pyrrolidine rings, each closed via its
    # own CD-N bond cut into a loop-closure DistanceConstraint (World.cpp
    # ~:679; spec Sec. 3.2). SHALL: >= 1; the structural fact for THIS
    # molecule is exactly 2 -- assert both, so a silent drop to 0 or 1 is
    # caught.
    assert world.num_loop_constraints >= 1, (
        "Ace-Pro-Pro-Nme must yield at least one ring-closing loop constraint "
        "-- otherwise the cyclic Fixman path is untested/vacuous"
    )
    assert world.num_loop_constraints == 2, (
        f"Ace-Pro-Pro-Nme is DI-proline (2 pyrrolidine rings): expected 2 loop "
        f"constraints, got {world.num_loop_constraints}"
    )

    ld_start = world.current_constraint_log_det()
    assert ld_start != 0.0, (
        "calcConstraintLogDet is 0 at the starting geometry -- the loop term "
        "never fired (T2.0 would be vacuous)"
    )
    assert np.isfinite(ld_start)

    # Conformation dependence: move to a different geometry (a short,
    # AlwaysAccept-equivalent MD-HMC burn-in under Fixman) and confirm the
    # loop-closure log-det actually CHANGES -- catches a calcConstraintLogDet
    # that returns a config-INDEPENDENT constant (which would make the SHALL
    # above trivially true but the term itself useless).
    ctx.run_rex(0, 40, 0, False)
    ld_after = world.current_constraint_log_det()
    assert np.isfinite(ld_after)
    assert abs(ld_after - ld_start) > 1e-6, (
        f"calcConstraintLogDet did not change across 40 MD-HMC rounds "
        f"({ld_start} -> {ld_after}) -- expected conformation-dependence"
    )


# ---------------------------------------------------------------------------
#  T2.1 -- reference recovery (slow). SHALL (ordering) + SHOULD (absolute
#  threshold, calibrated from this run and recorded below).
# ---------------------------------------------------------------------------
def test_ace_pro_pro_nme_fixman_closer_to_flexible_than_torsional(tmp_path):
    if not slow_tests_enabled():
        pytest.skip("slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)")
    _skip_if_missing()

    n_equil, n_prod, write_freq, nbins = 300, 1500, 1, 12

    # Build once (bonds/quads are geometry-derived, independent of use_fixman).
    ctx_probe, _world_probe, bonds = build_torsional_context(
        str(tmp_path / "pro_probe"), PRMTOP, RST7, seed=200, use_fixman=True
    )
    quads = _dihedral_quads_prmtop_order(ctx_probe, bonds)
    assert len(quads) >= 5, "expected several backbone/chi torsions with valid quads"

    ctx_flex = build_cartesian_context(str(tmp_path / "pro_flex"), PRMTOP, RST7, seed=201)
    dih_flex = run_and_load_dihedrals(
        ctx_flex, str(tmp_path / "pro_flex"), PRMTOP, n_equil, n_prod, write_freq, quads
    )

    ctx_tor, _world_tor, _ = build_torsional_context(
        str(tmp_path / "pro_tor"), PRMTOP, RST7, seed=202, use_fixman=False, bonds=bonds
    )
    dih_tor = run_and_load_dihedrals(
        ctx_tor, str(tmp_path / "pro_tor"), PRMTOP, n_equil, n_prod, write_freq, quads
    )

    ctx_fix, _world_fix, _ = build_torsional_context(
        str(tmp_path / "pro_fix"), PRMTOP, RST7, seed=203, use_fixman=True, bonds=bonds
    )
    dih_fix = run_and_load_dihedrals(
        ctx_fix, str(tmp_path / "pro_fix"), PRMTOP, n_equil, n_prod, write_freq, quads
    )

    # Hellinger on IACT-THINNED series (each world/torsion thinned by its own
    # integrated autocorrelation time). Robosample MD-HMC output is heavily
    # autocorrelated; comparing raw per-round frames measures noise, not the
    # ensembles (this was the original T2.1 failure -- see the idealized-chain
    # module docstring). A torsion is INFORMATIVE only if its FLEXIBLE reference
    # actually explored (>= 4 populated bins); a near-frozen torsion carries no
    # signal about the FIX-vs-TOR ordering and is excluded from the comparison.
    def _thinned_counts(series):
        thin, _tau, _s = thin_to_independent(series)
        c, _ = histogram_counts(thin, nbins)
        return c.astype(float)

    h_fix_flex = []
    h_tor_flex = []
    for k in range(len(quads)):
        c_flex = _thinned_counts(dih_flex[:, k])
        if int(np.count_nonzero(c_flex)) < 4:
            print(f"[pro] torsion {k}: FLEX near-frozen ({int(np.count_nonzero(c_flex))} bins) -- excluded")
            continue
        c_tor = _thinned_counts(dih_tor[:, k])
        c_fix = _thinned_counts(dih_fix[:, k])
        hf = hellinger_distance(c_fix, c_flex)
        ht = hellinger_distance(c_tor, c_flex)
        h_fix_flex.append(hf)
        h_tor_flex.append(ht)
        print(f"[pro] torsion {k}: H(FIX,FLEX)={hf:.4f}  H(TOR,FLEX)={ht:.4f}")

    assert len(h_fix_flex) >= 3, (
        f"only {len(h_fix_flex)} informative torsions -- proline sampling too poor to judge the ordering"
    )
    h_fix_flex = np.asarray(h_fix_flex)
    h_tor_flex = np.asarray(h_tor_flex)
    mean_h_fix_flex = float(np.mean(h_fix_flex))
    mean_h_tor_flex = float(np.mean(h_tor_flex))
    closer = int(np.sum(h_fix_flex < h_tor_flex))
    print(
        f"[pro] mean H(FIX,FLEX)={mean_h_fix_flex:.4f}  mean H(TOR,FLEX)={mean_h_tor_flex:.4f}  "
        f"FIXMAN closer on {closer}/{len(h_fix_flex)}"
    )

    # SHALL: FIXMAN closer to FLEXIBLE than TORSIONAL (Kandel 2016 qualitative
    # ordering). Judged on the MAJORITY of informative torsions AND in the mean,
    # so a single noisy torsion cannot flip the verdict and vacuum's residual
    # extrinsic 1-4 coupling (which Fixman provably cannot remove, spec Sec. 3.2
    # NOTE) does not force a false failure on any one angle.
    assert closer > len(h_fix_flex) / 2 and mean_h_fix_flex < mean_h_tor_flex, (
        f"FIXMAN should be closer to FLEXIBLE than TORSIONAL: closer on "
        f"{closer}/{len(h_fix_flex)}, mean H(FIX,FLEX)={mean_h_fix_flex:.4f} "
        f"vs H(TOR,FLEX)={mean_h_tor_flex:.4f}"
    )
