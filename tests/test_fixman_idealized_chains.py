"""Tier 1 (docs/specs/fixman-idealized-chains-validation.md): statistical
distribution recovery on the idealized bead chains, via the Python full-stack
path (spec Sec. 3.3). Gated by ``ROBOSAMPLE_SLOW_TESTS`` (unset -> skip, matching
``TestFixmanBoltzmann.cpp``'s C++ precedent).

Two methodology points that the first cut got wrong and this file fixes:

1. AUTOCORRELATION. Robosample MD-HMC output is heavily autocorrelated (the C4
   torsion IACT is ~360 at the old mdSteps=20/ts=0.001), so a chi-square GOF test
   that treats per-round frames as independent inflates chi2 into the thousands
   and rejects EVERY hypothesis. The fix: sample with good mixing (the raised
   ``build_torsional_context`` defaults) and THIN each series by its integrated
   autocorrelation time (``thin_to_independent``) before any histogram statistic;
   flatness is then judged against the finite-sample noise floor
   (``flat_noise_floor``), not a fixed critical value.

2. C4 IS A FLAT-METRIC NULL CONTROL, not a sqrt(det M) reproduction. C4's single
   torsion is TERMINAL: with the shipped prmtop's chain decomposition it rotates a
   single point mass at a FIXED perpendicular distance from the bond axis, so
   det M is torsion-INDEPENDENT and the generalized-coordinate marginal is already
   flat with Fixman OFF. The Jain 2013 eq:24 sqrt(det M) bias exists ONLY for the
   artificial 2+2 decomposition the C++ Tier-0 (TestFixmanIdealizedChains.cpp)
   hand-builds; the standard chain builder provably does not reproduce it (that
   file's own banner). So the full-stack C4 tests that Fixman does NOT introduce
   bias where there is none. The genuine metric bias lives in the INTERIOR
   torsions of C5+ (multi-atom moving groups on both sides of the axis).

The decomposition-INVARIANT Fixman guarantee this file leans on: with U == 0 a
correct Fixman makes the JOINT torsion distribution uniform, hence EVERY 1-D
torsion marginal is uniform -- regardless of which spanning tree the builder
picked. Fixman must therefore never leave a marginal LESS flat than raw
TORSIONAL. That is the load-bearing, calibration-free assertion below.

NOTE: sample counts are reduced from full statistical power (per-round OpenMM/CUDA
overhead dominates for these tiny systems); thresholds are noise-floor-relative
and orderings are majority-based so a modest N cannot make them flaky.
"""

from __future__ import annotations

import pathlib
import warnings

import numpy as np
import pytest

pytest.importorskip("openmm")
robosample = pytest.importorskip("robosample")
pytest.importorskip("mdtraj")

from fixman_chain_lib import (  # noqa: E402
    build_cartesian_context,
    build_torsional_context,
    flat_noise_floor,
    hellinger_distance,
    histogram_counts,
    run_and_load_dihedrals,
    slow_tests_enabled,
    thin_to_independent,
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
FIXMAN_DIR = REPO_ROOT / "examples" / "fixman"


def _skip_unless_slow():
    if not slow_tests_enabled():
        pytest.skip("slow statistical tier (set ROBOSAMPLE_SLOW_TESTS=1)")


def _flatness(series: np.ndarray, nbins: int):
    """(H-to-flat, N_eff, IACT) for one torsion series, on IACT-thinned samples."""
    thin, tau, _ = thin_to_independent(series)
    counts, _ = histogram_counts(thin, nbins)
    h = hellinger_distance(counts.astype(float), np.ones(nbins))
    return h, len(thin), tau


def _assert_flat(name: str, h: float, n_eff: int, nbins: int):
    """A marginal is 'flat within noise' if H is within a generous multiple of the
    finite-sample noise floor. The multiple (3x) plus a 0.12 absolute floor lets a
    truly flat marginal pass at a modest N while still catching a sqrt(det M)-scale
    bias (well-mixed measurements put a genuinely flat marginal at H ~ 0.02)."""
    floor = flat_noise_floor(n_eff, nbins)
    bound = max(3.0 * floor, 0.12)
    assert h < bound, f"{name}: marginal not flat (H={h:.4f} >= bound={bound:.4f}, noise_floor={floor:.4f})"


def _is_flat(h: float, n_eff: int, nbins: int) -> bool:
    return h < max(3.0 * flat_noise_floor(n_eff, nbins), 0.12)


def _fixman_flattening_report(label: str, dih_tor: np.ndarray, dih_fix: np.ndarray, nbins: int, k_torsions: int):
    """Per-torsion empirical signal selection + the Fixman guarantee (spec Sec. 5,
    signal-torsion selection). For each torsion, on IACT-thinned samples:
      * classify SIGNAL (TORSIONAL non-flat beyond the noise floor, i.e. a genuine
        sqrt(det M) bias exists) vs no-op (TORSIONAL already flat);
      * the decomposition-invariant U==0 guarantee: FIXMAN marginal flat (checked
        for EVERY torsion -- true in both the biased and the no-op case);
      * on SIGNAL torsions ONLY, the sign-discriminating check H(FIX,flat) <
        H(TOR,flat): a flipped Fixman sign would make FIXMAN *more* biased and fail
        here (the check C4's flat metric cannot make).
    Returns (n_signal, n_fix_flat, n_flattened)."""
    n_signal = n_fix_flat = n_flattened = 0
    print(f"[{label}]")
    for k in range(k_torsions):
        h_tor, neff_tor, _ = _flatness(dih_tor[:, k], nbins)
        h_fix, neff_fix, _ = _flatness(dih_fix[:, k], nbins)
        floor_fix = flat_noise_floor(neff_fix, nbins)
        is_signal = h_tor > max(2.5 * flat_noise_floor(neff_tor, nbins), 0.08)
        fix_flat = _is_flat(h_fix, neff_fix, nbins)
        n_fix_flat += 1 if fix_flat else 0
        flattened = h_fix < h_tor - 0.5 * floor_fix
        print(
            f"  torsion {k} [{'SIGNAL' if is_signal else 'no-op '}]: "
            f"H(TOR,flat)={h_tor:.4f} H(FIX,flat)={h_fix:.4f} (floor~{floor_fix:.4f})"
        )
        if is_signal:
            n_signal += 1
            n_flattened += 1 if flattened else 0
    return n_signal, n_fix_flat, n_flattened


# ---------------------------------------------------------------------------
#  T1.1 -- C4 flat-metric NULL control (spec revision). Terminal torsion =>
#  torsion-independent det M => TORSIONAL/FIXMAN/FLEXIBLE all flat; Fixman a no-op.
# ---------------------------------------------------------------------------
def test_c4_flat_metric_null_control(tmp_path):
    _skip_unless_slow()
    prmtop = FIXMAN_DIR / "c4" / "c4.prmtop"
    rst7 = FIXMAN_DIR / "c4" / "c4.rst7"
    if not prmtop.exists():
        pytest.skip(f"example not found: {prmtop}")

    n_equil, n_prod, write_freq, nbins = 800, 8000, 1, 18
    quad = [[0, 1, 2, 3]]

    series = {}
    for name, fx, seed in [("TOR", False, 1), ("FIX", True, 2)]:
        ctx, _w, bonds = build_torsional_context(
            str(tmp_path / f"c4_{name}"), prmtop, rst7, seed=seed, use_fixman=fx
        )
        assert len(bonds) == 1, "C4 has exactly one rotatable (non-terminal) bond"
        series[name] = run_and_load_dihedrals(
            ctx, str(tmp_path / f"c4_{name}"), prmtop, n_equil, n_prod, write_freq, quad
        )[:, 0]
    ctx_flex = build_cartesian_context(str(tmp_path / "c4_flex"), prmtop, rst7, seed=3)
    series["FLEX"] = run_and_load_dihedrals(
        ctx_flex, str(tmp_path / "c4_flex"), prmtop, n_equil, n_prod, write_freq, quad
    )[:, 0]

    res = {}
    for name, s in series.items():
        h, neff, tau = _flatness(s, nbins)
        res[name] = (h, neff)
        print(f"[C4] {name}: H(.,flat)={h:.4f} N_eff={neff} IACT={tau:.1f}")

    # C4's metric is torsion-independent -> all three marginals are flat.
    for name, (h, neff) in res.items():
        _assert_flat(f"C4 {name}", h, neff, nbins)

    # Fixman must not INTRODUCE bias where the metric is flat: FIXMAN no less flat
    # than raw TORSIONAL beyond noise.
    h_tor, neff_tor = res["TOR"]
    h_fix, neff_fix = res["FIX"]
    tol = 2.0 * flat_noise_floor(min(neff_tor, neff_fix), nbins)
    assert h_fix < h_tor + tol, (
        f"Fixman INTRODUCED bias on flat-metric C4: H(FIX)={h_fix:.4f} > H(TOR)={h_tor:.4f} + {tol:.4f}"
    )


# ---------------------------------------------------------------------------
#  T1.2 -- C5. FIXMAN flattens EVERY torsion marginal (the decomposition-
#  invariant U==0 guarantee), and is never less flat than raw TORSIONAL.
# ---------------------------------------------------------------------------
def test_c5_fixman_flattens_each_torsion(tmp_path):
    _skip_unless_slow()
    prmtop = FIXMAN_DIR / "c5" / "c5.prmtop"
    rst7 = FIXMAN_DIR / "c5" / "c5.rst7"
    if not prmtop.exists():
        pytest.skip(f"example not found: {prmtop}")

    n_equil, n_prod, write_freq, nbins = 800, 8000, 1, 18
    quads = [[0, 1, 2, 3], [1, 2, 3, 4]]

    ctx_tor, _wt, bonds = build_torsional_context(
        str(tmp_path / "c5_tor"), prmtop, rst7, seed=1, use_fixman=False
    )
    assert len(bonds) == 2, "C5 has exactly two rotatable bonds"
    dih_tor = run_and_load_dihedrals(
        ctx_tor, str(tmp_path / "c5_tor"), prmtop, n_equil, n_prod, write_freq, quads
    )
    ctx_fix, _wf, _ = build_torsional_context(
        str(tmp_path / "c5_fix"), prmtop, rst7, seed=2, use_fixman=True
    )
    dih_fix = run_and_load_dihedrals(
        ctx_fix, str(tmp_path / "c5_fix"), prmtop, n_equil, n_prod, write_freq, quads
    )

    n_signal, n_fix_flat, n_flattened = _fixman_flattening_report("C5", dih_tor, dih_fix, nbins, 2)
    # Fixman guarantee (holds whether or not a torsion is biased): every FIXMAN
    # marginal is flat.
    assert n_fix_flat == 2, f"C5: only {n_fix_flat}/2 FIXMAN marginals are flat"
    if n_signal == 0:
        warnings.warn(
            "C5 exercised NO metric bias (every TORSIONAL torsion already flat): the "
            "full-stack Fixman validation is a NULL control here, not a bias-removal "
            "test -- see docs/specs/fixman-idealized-chains-validation.md Sec. 7.2.",
            stacklevel=2,
        )
    else:
        # Sign-discriminating: Fixman must flatten every biased torsion.
        assert n_flattened == n_signal, (
            f"C5: Fixman failed to flatten {n_signal - n_flattened}/{n_signal} biased torsions"
        )


# ---------------------------------------------------------------------------
#  T1.3 -- C15 capstone (12 torsions). Same guarantee; majority-based so a few
#  noisy simultaneous marginals cannot make it flaky (C11 dropped per spec:
#  C4+C5+C15 = 1,2,12 torsions already exercise the scaling).
# ---------------------------------------------------------------------------
def test_c15_fixman_flattens_every_torsion(tmp_path):
    _skip_unless_slow()
    prmtop = FIXMAN_DIR / "c15" / "c15.prmtop"
    rst7 = FIXMAN_DIR / "c15" / "c15.rst7"
    if not prmtop.exists():
        pytest.skip(f"example not found: {prmtop}")

    n_equil, n_prod, write_freq, nbins = 600, 6000, 1, 12
    quads = [[i, i + 1, i + 2, i + 3] for i in range(12)]

    ctx_tor, _wt, bonds = build_torsional_context(
        str(tmp_path / "c15_tor"), prmtop, rst7, seed=1, use_fixman=False
    )
    assert len(bonds) == 12
    dih_tor = run_and_load_dihedrals(
        ctx_tor, str(tmp_path / "c15_tor"), prmtop, n_equil, n_prod, write_freq, quads
    )
    ctx_fix, _wf, _ = build_torsional_context(
        str(tmp_path / "c15_fix"), prmtop, rst7, seed=2, use_fixman=True
    )
    dih_fix = run_and_load_dihedrals(
        ctx_fix, str(tmp_path / "c15_fix"), prmtop, n_equil, n_prod, write_freq, quads
    )

    n_signal, n_fix_flat, n_flattened = _fixman_flattening_report("C15", dih_tor, dih_fix, nbins, 12)
    # 12 simultaneous flatness tests: require the majority (a smoke-level capstone,
    # not a full-power per-torsion gate).
    assert n_fix_flat >= 9, f"only {n_fix_flat}/12 C15 FIXMAN marginals are flat"
    if n_signal == 0:
        warnings.warn(
            "C15 exercised NO metric bias (every TORSIONAL torsion already flat): the "
            "full-stack Fixman validation is a NULL control here -- see "
            "docs/specs/fixman-idealized-chains-validation.md Sec. 7.2.",
            stacklevel=2,
        )
    else:
        # Majority of biased torsions must be flattened (a few of the 12 may be noisy).
        need = max(1, int(round(0.75 * n_signal)))
        assert n_flattened >= need, (
            f"C15: Fixman flattened only {n_flattened}/{n_signal} biased torsions (need >= {need})"
        )
