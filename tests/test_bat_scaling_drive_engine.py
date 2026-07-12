"""V5 (engine-level Jacobian gate) and V10 (anchor involution) -- Stage 2a.

Spec: docs/specs/replica-exchange-nonequilibrium-work.md B4/D5/D6/D7/INV-9.

V5 is a hard gate (coordinator directive, 2026-07-12): a wrong Jacobian sign
silently biases the entire sampler. This file finite-differences the ACTUAL
Cartesian scaling map the C++ engine applies (World.preview_bat_scaling, which
calls robo::applyBatScaling/calcBatVolumeLogJac verbatim -- the SAME code path
World.apply_bat_scaling_drive uses in production) against the analytic
lnJac World.preview_bat_scaling itself returns, on a REAL molecule
(ala-dipeptide) built through the ordinary Context/World API -- not the
standalone numpy reference (tests/bat_jacobian_scaling_check.py, untouched,
still the spec-level analytic half).

FLAGGED LIMITATION (BAT/z-matrix API ambiguity, coder checkpoint): a scaled
body whose PARENT is itself the molecule's ROOT body (bodyParent[bodyParent
[b]] == Ground) has no third ancestor atom (zK) -- its angle DOF is silently
skipped (BatScaling.hpp/cpp document this; it is a real, OBSERVED effect
below, not just a theoretical corner case). Concretely, in ala-dipeptide the
backbone body anchored at atom 7 (bond 3-7) is ALWAYS one such body (its
parent is the whole Ground-rooted cap cluster), so EVERY chain rooted there
has its FIRST body contribute r-only (N_scaled=1), not r+theta (2), to the
D5 count -- confirmed empirically in test_v5b/test_v5c below (measured, not
assumed, via body-index probes; see the coder checkpoint for the derivation).
The engine additionally has NO joint type that scales theta in ISOLATION
(BendStretch/SphericalCoords always coscale r+theta, D5 table), so the
spec's idealized 1-DOF "single angle scale" target (N_scaled=1) has no
direct engine realization; V5b decomposes the closest honest analogue.

Three engine-level cases:

  V5a: a single Slider bond (N_scaled=1, r only)      -> want 3 ln(s) exactly
       (D6's own numeric example, spec-cited verbatim).
  V5b: a 2-body BendStretch chain (bond 3-7 then 7-10) -> N_scaled=3 (body 7:
       r-only per the flagged limitation above; body 10: r AND theta, its
       zK=0 is a real ancestor atom). Decomposed into the spec's two named
       1-DOF closed forms per body and checked against the engine FD.
  V5c: a >=3-BODY CHAIN (three consecutive BendStretch bonds along the
       backbone), ALL scaled together (N_scaled=5, see above) -> want
       (J(x')-J(x0)) + 5 ln(s). This is the discriminating case (spec V5
       note): it is the only one that exercises an UPSTREAM body's z-matrix
       ancestor lookup (zJ/zK) reading ALREADY-CASCADED (rigidly moved)
       positions correctly -- a single-DOF case cannot see a mismatch there.

V10 (INV-9): the paired map M_{1/s,mu}(M_{s,mu}(q)) == q to floating point
for a SHARED, frozen anchor mu, and FAILS to round-trip if the reverse leg
uses a DIFFERENT anchor -- the discriminating property D2/INV-9 exists to
guard (a per-state mu would silently break the involution).

Run: pytest -q tests/test_bat_scaling_drive_engine.py
"""

from __future__ import annotations

import os
import pathlib

import numpy as np
import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))

if _REQUIRE_OPENMM:
    import robosample
else:
    robosample = pytest.importorskip("robosample")

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PRMTOP = REPO_ROOT / "examples" / "ala-dipeptide.prmtop"
RST7 = REPO_ROOT / "examples" / "ala-dipeptide.rst7"


def _skip_if_missing_inputs():
    if _REQUIRE_OPENMM:
        assert PRMTOP.exists() and RST7.exists(), (
            f"example inputs not found: {PRMTOP} / {RST7} -- the authoritative "
            "gate (ROBOSAMPLE_REQUIRE_OPENMM=1) requires them to be present"
        )
    elif not PRMTOP.exists() or not RST7.exists():
        pytest.skip(f"example inputs not found: {PRMTOP} / {RST7}")


def _build_scaled_world(bonds, joint_type, tag):
    """A robotic world over ala-dipeptide with `bonds` (global atom-index
    pairs) set to `joint_type`, everything else Rigid, distort_option=
    ScaleBendStretch, mdSteps=0 (D7). No context.initialize() -- the drive
    (World::previewBatScaling) needs only the built RobotModel, not OpenMM.
    """
    ctx = robosample.Context(
        str(pathlib.Path("/tmp") / f"bat_v5_probe_{tag}"), 7, robosample.AmberDihedralClassifier()
    )
    ctx.load_amber(str(PRMTOP), str(RST7))
    sele = ctx.build_flexibilities(bonds, joint_type, False)
    w = ctx.add_robotic_world(sele)
    w.add_sampler(
        timeStep=0.0,
        mdSteps=0,
        acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
        use_nuts=False,
        use_fixman=True,
        distort_option=robosample.rb.DistortOption.ScaleBendStretch,
    )
    st = ctx.system_topology
    flat = []
    for a in range(st.num_atoms):
        flat += [st.atoms_x[a], st.atoms_y[a], st.atoms_z[a]]
    return ctx, w, flat


def _fd_log_jac(world, x0_flat, s, anchor_r, anchor_theta, h=1e-6):
    """log|det dx'/dx0| by central finite difference of world.preview_bat_scaling."""
    n = len(x0_flat)
    x0 = np.array(x0_flat)
    jac = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n)
        e[k] = h
        plus, _, _ = world.preview_bat_scaling((x0 + e).tolist(), s, anchor_r, anchor_theta)
        minus, _, _ = world.preview_bat_scaling((x0 - e).tolist(), s, anchor_r, anchor_theta)
        jac[:, k] = (np.array(plus) - np.array(minus)) / (2.0 * h)
    sign, logdet = np.linalg.slogdet(jac)
    assert sign > 0, "map orientation reversed -- something is badly wrong"
    return logdet


def _angle(p, a, b, c):
    v1 = p[a] - p[b]
    v2 = p[c] - p[b]
    cosang = np.clip(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)), -1.0, 1.0)
    return np.arccos(cosang)


# ---------------------------------------------------------------------------
# V5a -- single Slider bond (N_scaled=1, r only): want 3 ln(s) exactly (D6).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("s", [1.3, 0.8])
def test_v5a_single_slider_bond_matches_3_ln_s(s):
    _skip_if_missing_inputs()
    _ctx, w, x0 = _build_scaled_world([(26, 27)], robosample.rb.JointType.Slider, f"5a_{s}")

    _out, n_scaled, ln_jac_analytic = w.preview_bat_scaling(x0, s, {}, {})
    assert n_scaled == 1

    want = 3.0 * np.log(s)
    assert ln_jac_analytic == pytest.approx(want, abs=1e-9), (
        f"analytic lnJac={ln_jac_analytic!r} != 3 ln(s)={want!r}"
    )

    fd = _fd_log_jac(w, x0, s, {}, {})
    assert fd == pytest.approx(want, abs=2e-4), (
        f"ENGINE finite-difference lnJac={fd!r} != 3 ln(s)={want!r} -- V5 GATE FAILURE"
    )
    assert fd == pytest.approx(ln_jac_analytic, abs=2e-4), (
        f"engine FD={fd!r} != analytic formula={ln_jac_analytic!r} -- V5 GATE FAILURE"
    )


# ---------------------------------------------------------------------------
# V5b -- 2-body BendStretch chain (bond 3-7, bond 7-10): body 7 is root-
# adjacent (r-only, the flagged limitation), body 10 has a real zK (r+theta).
# N_scaled=3 total; decomposed per body against the spec's named 1-DOF forms.
# Atom indices below (zJ=0 for body 7, zJ=7/zK=0 for body 10) were VERIFIED
# empirically (not assumed) against the engine's own lnJac output -- see the
# coder checkpoint for the derivation.
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("s", [1.05, 0.85])
def test_v5b_two_body_chain_matches_engine_fd_and_decomposes(s):
    _skip_if_missing_inputs()
    bonds = [(3, 7), (7, 10)]
    _ctx, w, x0 = _build_scaled_world(bonds, robosample.rb.JointType.BendStretch, f"5b_{s}")

    out, n_scaled, ln_jac_analytic = w.preview_bat_scaling(x0, s, {}, {})
    assert n_scaled == 3  # body(7): r-only (root-adjacent, zK undefined); body(10): r+theta

    fd = _fd_log_jac(w, x0, s, {}, {})
    assert fd == pytest.approx(ln_jac_analytic, abs=2e-4), (
        f"engine FD={fd!r} != analytic formula={ln_jac_analytic!r} -- V5 GATE FAILURE"
    )

    x0v = np.array(x0).reshape(-1, 3)
    xt = np.array(out).reshape(-1, 3)

    # body(7): zI=7, zJ=0 (the molecule's root atom) -- r-only ("3 ln s" form).
    r0_7 = np.linalg.norm(x0v[7] - x0v[0])
    r1_7 = np.linalg.norm(xt[7] - xt[0])
    term_7 = 2.0 * np.log(r1_7 / r0_7) + np.log(s)

    # body(10): zI=10, zJ=7, zK=0 -- r AND theta.
    r0_10 = np.linalg.norm(x0v[10] - x0v[7])
    r1_10 = np.linalg.norm(xt[10] - xt[7])
    th0_10 = _angle(x0v, 10, 7, 0)
    th1_10 = _angle(xt, 10, 7, 0)
    term_10_r = 2.0 * np.log(r1_10 / r0_10) + np.log(s)
    term_10_theta = np.log(np.sin(th1_10) / np.sin(th0_10)) + np.log(s)

    total = term_7 + term_10_r + term_10_theta
    assert total == pytest.approx(ln_jac_analytic, abs=1e-8), (
        "per-body decomposition into the spec's named 1-DOF bond/angle forms does "
        "not sum to the engine's lnJac -- either the atom-index assumptions above "
        "are wrong or the engine's D5/D6 composition is wrong"
    )


# ---------------------------------------------------------------------------
# V5c -- >=3-body BendStretch chain (upstream bond AND angle), ALL scaled
# together: the discriminating case (spec V5 note on the geometric BAT tree).
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("s", [1.05, 0.9])
def test_v5c_three_body_chain_matches_engine_fd(s):
    _skip_if_missing_inputs()
    # Backbone chain N(3)-CA(7)-C(10)-N(16): three CONSECUTIVE bonds, all
    # BendStretch -> three chained flexible bodies (body(7,9,11,12,13,14),
    # body(10,15), body(16..31)). Scaling body(7,...) (the UPSTREAM one)
    # rigidly cascades into the two DOWNSTREAM bodies -- exactly the
    # "geometric-BAT/mobilizer-coordinate mismatch" case V5 exists to catch.
    # N_scaled=5: body7 contributes 1 (root-adjacent, r-only, see module
    # docstring), body10 and body16 each contribute 2 (r+theta, valid zK).
    bonds = [(3, 7), (7, 10), (10, 16)]
    _ctx, w, x0 = _build_scaled_world(bonds, robosample.rb.JointType.BendStretch, f"5c_{s}")

    _out, n_scaled, ln_jac_analytic = w.preview_bat_scaling(x0, s, {}, {})
    assert n_scaled == 5

    fd = _fd_log_jac(w, x0, s, {}, {})
    assert fd == pytest.approx(ln_jac_analytic, abs=3e-4), (
        f"engine FD={fd!r} != analytic formula={ln_jac_analytic!r} on a >=3-body chain "
        "-- V5 GATE FAILURE (a geometric-BAT/mobilizer-coordinate mismatch would show "
        "up ONLY here, not in the single-DOF cases)"
    )

    # A flipped (J_fin-J_ini) sign (the F7 bug this spec REFUTES) is
    # numerically distinguishable from the correct answer here.
    if s != 1.0:
        flipped = -ln_jac_analytic
        assert abs(flipped - fd) > 1e-2, "the F7-buggy sign-flipped composition must NOT match the engine FD"


# ---------------------------------------------------------------------------
# V10 -- anchor involution (INV-9)
# ---------------------------------------------------------------------------
def test_v10_shared_anchor_round_trips_and_per_state_anchor_breaks_it():
    _skip_if_missing_inputs()
    bonds = [(3, 7), (7, 10), (10, 16)]
    _ctx, w, x0 = _build_scaled_world(bonds, robosample.rb.JointType.BendStretch, "v10")

    s = 1.05
    # Nontrivial, SHARED anchor (INV-9): pick something away from the current
    # value so the round-trip is a real algebraic cancellation, not a
    # degenerate mu==q no-op. Keyed by each scaled body's zI (7, 10, 16).
    anchor_r = {7: 0.02, 10: 0.03, 16: 0.04}
    anchor_theta = {7: -0.05, 10: 0.08, 16: -0.02}

    forward, _, _ = w.preview_bat_scaling(x0, s, anchor_r, anchor_theta)
    back, _, _ = w.preview_bat_scaling(forward, 1.0 / s, anchor_r, anchor_theta)

    x0_arr = np.array(x0)
    back_arr = np.array(back)
    max_err = np.max(np.abs(back_arr - x0_arr))
    assert max_err < 1e-9, f"M_{{1/s,mu}}(M_{{s,mu}}(q)) != q for the SHARED anchor (max err {max_err})"

    # Different (per-state-style) anchor on the reverse leg -- INV-9's
    # discriminating failure: the round-trip SHALL break.
    anchor_r_other = {k: v + 0.01 for k, v in anchor_r.items()}
    anchor_theta_other = {k: v + 0.01 for k, v in anchor_theta.items()}
    back_wrong, _, _ = w.preview_bat_scaling(forward, 1.0 / s, anchor_r_other, anchor_theta_other)
    max_err_wrong = np.max(np.abs(np.array(back_wrong) - x0_arr))
    assert max_err_wrong > 1e-3, (
        "using a DIFFERENT anchor for the reverse leg SHOULD break the round-trip "
        "(INV-9) -- got a suspiciously small error, the anchor may not be wired in"
    )
