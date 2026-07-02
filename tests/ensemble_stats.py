"""ensemble_stats.py -- shared statistics for the ensemble-validation suite.

Implements the estimators and hygiene rules from
``docs/specs/ensemble-validation/00-foundations.md`` (see especially
Sec. 5.1 "Shirts ensemble-validation estimators" and Sec. 6 "Statistical
hygiene"). Both Tier 1 (``test_ensemble_pe_ladder.py``) and Tier 2
(``test_torsion_conformational.py``) import this module; nothing here is
Robosample-specific beyond the ``moves.csv`` column layout, so it is pure
numpy/scipy/pandas.

This module is ALSO a pytest file: the ``test_*`` functions at the bottom are
the mandated self-check that validates the STATISTIC (not the sampler)
against synthetic data with a known, closed-form answer. They are always-on
(no ``ROBOSAMPLE_SLOW_TESTS`` gate) because they run in well under a second
and are what guards the estimator code itself; ``python3 -m pytest
tests/ensemble_stats.py -q`` runs only these.

Conventions fixed here (foundations Sec. 5.1 "sign discipline")
-----------------------------------------------------------------
* ``y_k = ln(n_{k,2} / n_{k,1})`` (dataset 2 over dataset 1) is the response
  in Estimator A; the expected slope is ``-(beta2 - beta1)`` where
  ``beta_i = 1/(kB*T_i)``. Estimator B uses the same sign: label 1 (=dataset
  2, i.e. logit increases with dataset-2 membership) is the positive class.
  This is the ``eq:6`` convention (``P2/P1``, not ``eq:7``'s ``P1/P2``); every
  caller in this suite MUST pass ``(e1, T1, e2, T2)`` in that fixed order.
* ``kB`` is in kJ/mol/K (matches ``tests/StatTest.hpp`` / ``TestEnsembleValidation.cpp``).

N_eff usage (foundations Sec. 6)
---------------------------------
Every error bar and every Shirts weight in this module derives from the
effective sample size ``N_eff = N/g``, ``g`` the statistical inefficiency
(Sec. 6), NOT the raw sample count ``N``. Two different, individually
standard techniques are used, chosen per statistic:

* **Mean/variance comparisons** (``mean_stderr_neff``): keep the full-N
  sample mean (unbiased, more precise) but divide its variance by ``N_eff``
  instead of ``N`` -- the standard block-averaging correction.
* **Histogram/MLE-based statistics** (KS test, Shirts A histogram counts,
  Shirts B logistic-regression Fisher information): systematically
  subsample (thin) the raw series with stride ``round(g)`` via
  :func:`thin_to_effective`, so the ``N`` that flows into the downstream
  formula (histogram counts, per-sample Fisher information) is honestly
  ``~N_eff`` roughly-independent draws, instead of algebraically rescaling a
  correlated-sample formula. Thinning is simpler to reason about and to
  test, at some cost in efficiency (some information is discarded) -- an
  explicit, surgical choice documented here rather than silently assumed.

``python/robosample/batstat.py`` and ``autoblock.py`` (named in the task as
candidate N_eff providers) turn out to be MDAnalysis/community-detection
dihedral-correlation tooling with no scalar-series statistical-inefficiency
estimator to reuse; ``statistical_inefficiency`` below implements the
standard Geyer (1992) initial-positive-sequence (IPS) estimator directly
(self-checked against the closed-form AR(1) result below) rather than
force-fitting an unrelated module.
"""

from __future__ import annotations

import dataclasses
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

# kB in kJ/mol/K (matches tests/StatTest.hpp / TestEnsembleValidation.cpp: the
# molar gas constant R, since energies here are per-mole).
KB_KJ_PER_MOL_K = 0.0083144626

MOVES_CSV_COLUMNS = [
    "round",
    "replica",
    "T",
    "world",
    "type",
    "kick",
    "accepted",
    "PE",
    "KE",
    "fixman",
    "H",
]


# ---------------------------------------------------------------------------
# moves.csv parsing (Context::runREX / Context.cpp, header fixed there)
# ---------------------------------------------------------------------------


def parse_moves_csv(path) -> pd.DataFrame:
    """Parse a ``{base}.moves.csv`` file into a DataFrame.

    Raises ``ValueError`` if any of the fixed columns
    (``round,replica,T,world,type,kick,accepted,PE,KE,fixman,H``) is missing
    -- a header mismatch means the C++ writer (``Context.cpp``) changed and
    this parser is now silently misreading columns, which must fail loudly
    rather than return garbage.
    """
    df = pd.read_csv(path)
    missing = set(MOVES_CSV_COLUMNS) - set(df.columns)
    if missing:
        raise ValueError(f"moves.csv {path} missing columns: {sorted(missing)}")
    return df


def world_series(
    df: pd.DataFrame,
    column: str,
    *,
    world_type: str | None = None,
    world_index: int | None = None,
    replica: int = 0,
    min_round: int = 0,
) -> np.ndarray:
    """Extract one world/replica's ``column`` (PE/KE/H/...) time series, in round order.

    ``min_round`` excludes equilibration rows: ``moves.csv`` logs every
    logged round regardless of equil/production phase (``Context.cpp``'s
    ``logThisRound`` gate is unconditional on the phase), so a caller that
    ran ``equil_rounds > 0`` MUST pass ``min_round=equil_rounds`` to drop the
    non-canonical burn-in rows before any statistic below is applied.
    """
    mask = (df["replica"] == replica) & (df["round"] >= min_round)
    if world_type is not None:
        mask &= df["type"] == world_type
    if world_index is not None:
        mask &= df["world"] == world_index
    sub = df.loc[mask].sort_values("round")
    if sub.empty:
        raise ValueError(
            f"No moves.csv rows for world_type={world_type!r} "
            f"world_index={world_index!r} replica={replica} min_round={min_round}"
        )
    return sub[column].to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# Autocorrelation-aware effective sample size (foundations Sec. 6)
# ---------------------------------------------------------------------------


def _normalized_autocovariance(x: np.ndarray, max_lag: int) -> np.ndarray:
    """rho[t] = autocovariance(t)/variance for t=0..max_lag; rho[0] == 1."""
    x = np.asarray(x, dtype=float)
    n = x.size
    xm = x - x.mean()
    var = np.dot(xm, xm) / n
    if var <= 0.0:
        return np.zeros(max_lag + 1)
    c = np.empty(max_lag + 1)
    for t in range(max_lag + 1):
        c[t] = np.dot(xm[: n - t], xm[t:]) / n
    return c / var


def statistical_inefficiency(x: np.ndarray, max_lag: int | None = None) -> float:
    """Geyer (1992) initial-positive-sequence estimate of ``g = 1 + 2*tau_int``.

    ``tau_int = sum_{t=1}^inf rho(t)`` (no 1/2 offset -- the convention under
    which ``N_eff = N/g = N/(1+2*sum rho(t))`` is the standard effective
    sample size, e.g. Sokal / Gelman-BDA3 / pymbar). Pairing lags as
    ``Gamma_m = rho(2m) + rho(2m+1)`` (Geyer's IPS trick keeps the partial
    sums monotonic under the AR-type processes this is used for) and summing
    while ``Gamma_m >= 0`` gives
    ``g = 2 * sum_{m=0}^{m*} Gamma_m - 1``.

    Self-checked below (``test_statistical_inefficiency_matches_ar1_theory``)
    against the closed-form AR(1) result ``g = (1+phi)/(1-phi)``.

    Returns 1.0 (no correction; conservative only in the sense of making
    NO assumption about correlation, so this is a caller-visible floor, not
    a silent no-op) when there are too few samples (<8) to estimate any lag.
    """
    x = np.asarray(x, dtype=float)
    n = x.size
    if n < 8:
        return 1.0
    if max_lag is None:
        max_lag = min(n - 1, 2000)
    max_lag = min(max_lag, n - 1)
    rho = _normalized_autocovariance(x, max_lag)
    total = 0.0
    m = 0
    while True:
        t0, t1 = 2 * m, 2 * m + 1
        if t1 > max_lag:
            break
        gamma = rho[t0] + rho[t1]
        if gamma < 0.0:
            break
        total += gamma
        m += 1
    g = 2.0 * total - 1.0
    return max(1.0, g)


def effective_sample_size(x: np.ndarray) -> float:
    """N_eff = N / g (foundations Sec. 6)."""
    x = np.asarray(x, dtype=float)
    g = statistical_inefficiency(x)
    return x.size / g


def thin_to_effective(x: np.ndarray) -> np.ndarray:
    """Systematically subsample ``x`` with stride ``round(g)`` (see module docstring).

    Yields ``~N_eff`` approximately-independent draws for statistics (KS,
    Shirts A/B) that need honestly-sized, roughly-iid input rather than an
    algebraic variance correction.
    """
    x = np.asarray(x, dtype=float)
    g = statistical_inefficiency(x)
    stride = max(1, int(round(g)))
    return x[::stride]


def mean_stderr_neff(x: np.ndarray) -> tuple[float, float]:
    """(mean, stderr) using the full-N mean but N_eff in the variance denominator."""
    x = np.asarray(x, dtype=float)
    n_eff = max(effective_sample_size(x), 1.0)
    var = float(np.var(x, ddof=1)) if x.size > 1 else 0.0
    return float(np.mean(x)), math.sqrt(var / n_eff)


def means_within_tolerance(
    x1: np.ndarray, x2: np.ndarray, k: float = 4.0
) -> tuple[bool, float, float]:
    """(ok, |mean1-mean2|, k*combined_stderr) -- the Sec. 6 moment pre-filter."""
    m1, se1 = mean_stderr_neff(x1)
    m2, se2 = mean_stderr_neff(x2)
    se = math.sqrt(se1**2 + se2**2)
    diff = abs(m1 - m2)
    return diff <= k * se, diff, k * se


# ---------------------------------------------------------------------------
# KS histogram equality (foundations Sec. 4/6: DOF-matched pairs only)
# ---------------------------------------------------------------------------


def ks_test_neff(x1: np.ndarray, x2: np.ndarray) -> tuple[float, float, int, int]:
    """Two-sample KS test on N_eff-thinned (~independent) samples.

    Returns (statistic, p_value, n1_used, n2_used).
    """
    from scipy import stats

    t1 = thin_to_effective(x1)
    t2 = thin_to_effective(x2)
    result = stats.ks_2samp(t1, t2)
    return float(result.statistic), float(result.pvalue), t1.size, t2.size


# ---------------------------------------------------------------------------
# Shirts two-temperature PE-slope test (foundations Sec. 5.1)
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class ShirtsFit:
    """Fitted Shirts PE-slope: ``slope`` estimates ``-(beta2-beta1)``."""

    slope: float
    se: float
    expected_slope: float
    z: float
    n_used: int


def shirts_slope_estimator_a(
    e1: np.ndarray,
    T1: float,
    e2: np.ndarray,
    T2: float,
    *,
    nbins: int = 40,
    min_count: float = 5.0,
    kB: float = KB_KJ_PER_MOL_K,
) -> ShirtsFit:
    """Estimator A: weighted-linear fit of the binned log-ratio (eq:6, eq:30).

    ``y_k = ln(n_{k,2}/n_{k,1})`` vs bin centre ``E_k``, weights
    ``w_k = 1/var(ln r_k)``, ``var(ln r_k) = 1/n_{k,1}+1/n_{k,2}-1/N1-1/N2``.
    Histogram counts are built on the N_eff-thinned series (module docstring)
    so ``n_{k,i}``/``N_i`` already reflect the effective sample size.
    """
    t1 = thin_to_effective(np.asarray(e1, dtype=float))
    t2 = thin_to_effective(np.asarray(e2, dtype=float))
    lo = min(t1.min(), t2.min())
    hi = max(t1.max(), t2.max())
    if hi <= lo:
        raise ValueError("Shirts estimator A: degenerate energy range, cannot bin")

    c1, edges = np.histogram(t1, bins=nbins, range=(lo, hi))
    c2, _ = np.histogram(t2, bins=nbins, range=(lo, hi))
    N1, N2 = float(t1.size), float(t2.size)
    centers = 0.5 * (edges[:-1] + edges[1:])

    Sw = Sx = Sy = Sxx = Sxy = 0.0
    used = 0
    for k in range(nbins):
        n1k, n2k = float(c1[k]), float(c2[k])
        if n1k < min_count or n2k < min_count:
            continue
        var = 1.0 / n1k + 1.0 / n2k - 1.0 / N1 - 1.0 / N2
        if var <= 0.0:
            continue
        w = 1.0 / var
        y = math.log(n2k / n1k)
        E = centers[k]
        Sw += w
        Sx += w * E
        Sy += w * y
        Sxx += w * E * E
        Sxy += w * E * y
        used += 1
    if used < 2:
        raise ValueError(
            "Shirts estimator A: fewer than 2 usable bins (widen the T gap or "
            "increase samples -- see foundations Sec. 6 temperature-gap guidance)"
        )
    denom = Sw * Sxx - Sx * Sx
    if denom <= 0.0:
        raise ValueError("Shirts estimator A: degenerate design matrix (denom<=0)")
    slope = (Sw * Sxy - Sx * Sy) / denom
    se = math.sqrt(Sw / denom)

    beta1 = 1.0 / (kB * T1)
    beta2 = 1.0 / (kB * T2)
    expected = -(beta2 - beta1)
    z = (slope - expected) / se if se > 0.0 else math.inf
    return ShirtsFit(slope=slope, se=se, expected_slope=expected, z=z, n_used=used)


def shirts_slope_estimator_b(
    e1: np.ndarray,
    T1: float,
    e2: np.ndarray,
    T2: float,
    *,
    kB: float = KB_KJ_PER_MOL_K,
    max_iter: int = 200,
    tol: float = 1e-10,
) -> ShirtsFit:
    """Estimator B: histogram-free ML logistic regression (eq:8, eq:36).

    ``ln L = sum_{T1} ln f(-a0-a1*E) + sum_{T2} ln f(a0+a1*E)``,
    ``f(x)=1/(1+exp(-x))`` -- equivalently ordinary logistic regression with
    label 1 for dataset-2 (T2) samples and 0 for dataset-1 (T1) samples.
    Concave -> unique maximum, found by Newton-Raphson/IRLS (the Fisher
    information IS the Hessian for a canonical-link GLM). ``se`` is
    ``sqrt(-1/(d2 lnL/d a1^2))`` (eq:36), i.e. the (a1,a1) entry of the
    inverse Fisher information.

    Fit on N_eff-thinned (~independent) samples (module docstring) so the
    naive per-sample Fisher information is not inflated by autocorrelation.
    Energies are standardized before the Newton iteration purely for
    numerical conditioning (kJ/mol magnitudes are large relative to 1);
    the returned slope/se are converted back to physical (1/energy) units.
    """
    t1 = thin_to_effective(np.asarray(e1, dtype=float))
    t2 = thin_to_effective(np.asarray(e2, dtype=float))
    E = np.concatenate([t1, t2])
    y = np.concatenate([np.zeros(t1.size), np.ones(t2.size)])

    E0 = E.mean()
    Es = E.std()
    if Es <= 0.0:
        raise ValueError("Shirts estimator B: zero-variance combined energy series")
    Ec = (E - E0) / Es
    X = np.column_stack([np.ones_like(Ec), Ec])

    beta = np.zeros(2)
    for _ in range(max_iter):
        eta = X @ beta
        p = 1.0 / (1.0 + np.exp(-eta))
        w = p * (1.0 - p)
        grad = X.T @ (y - p)
        fisher = X.T @ (X * w[:, None])
        try:
            step = np.linalg.solve(fisher, grad)
        except np.linalg.LinAlgError as exc:
            raise RuntimeError(
                "Shirts estimator B: singular Fisher information (degenerate data)"
            ) from exc
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    else:
        raise RuntimeError("Shirts estimator B: logistic MLE did not converge")

    eta = X @ beta
    p = 1.0 / (1.0 + np.exp(-eta))
    w = p * (1.0 - p)
    fisher = X.T @ (X * w[:, None])
    cov = np.linalg.inv(fisher)

    alpha1 = beta[1] / Es
    se = math.sqrt(cov[1, 1]) / Es

    beta1 = 1.0 / (kB * T1)
    beta2 = 1.0 / (kB * T2)
    expected = -(beta2 - beta1)
    z = (alpha1 - expected) / se if se > 0.0 else math.inf
    return ShirtsFit(slope=alpha1, se=se, expected_slope=expected, z=z, n_used=E.size)


# ===========================================================================
#  Always-on self-checks (no ROBOSAMPLE_SLOW_TESTS gate): validate the
#  STATISTIC against synthetic data with a known, closed-form answer, so the
#  estimator code itself is guarded independently of any real sampler run.
# ===========================================================================


def test_statistical_inefficiency_matches_ar1_theory():
    """g(AR(1) with lag-1 correlation phi) = (1+phi)/(1-phi), in closed form.

    A wrong sign, a missing "-1", or a factor-of-2 slip in
    ``statistical_inefficiency`` would miss this badly (e.g. phi=0.8 gives
    g_true=9; a "g=1+2*sum" bug without the IPS pairing correction, or a
    dropped "-1", is off by roughly a factor of 2) -- this is the
    discriminating check, not a re-statement of the implementation.
    """
    rng = np.random.default_rng(12345)
    for phi in (0.0, 0.3, 0.5, 0.8, 0.9):
        n = 200_000
        eps = rng.normal(size=n)
        x = np.empty(n)
        x[0] = eps[0]
        for i in range(1, n):
            x[i] = phi * x[i - 1] + eps[i]
        g_hat = statistical_inefficiency(x)
        g_true = (1.0 + phi) / (1.0 - phi)
        rel_err = abs(g_hat - g_true) / g_true
        assert rel_err < 0.15, (
            f"phi={phi}: g_true={g_true:.3f} g_hat={g_hat:.3f} rel_err={rel_err:.3%}"
        )


def test_effective_sample_size_of_iid_series_is_near_n():
    """An iid series has g~1, so N_eff ~ N (no spurious deflation)."""
    rng = np.random.default_rng(7)
    x = rng.normal(size=50_000)
    n_eff = effective_sample_size(x)
    assert 0.85 * x.size < n_eff <= x.size


def _synthetic_boltzmann_pe(rng, beta, n):
    """iid PE draws from density(E) ~ exp(-beta*E), E>=0 (flat density of states).

    The Shirts slope test is explicitly density-of-states-agnostic (foundations
    Sec. 5.1); an exponential is the simplest distribution with EXACTLY this
    Boltzmann form, so the true slope -(beta2-beta1) is known in closed form
    with no approximation, independent of the estimator under test.
    """
    return rng.exponential(scale=1.0 / beta, size=n)


def test_shirts_estimator_a_recovers_known_slope_on_synthetic_boltzmann_samples():
    rng = np.random.default_rng(2024)
    T1, T2 = 280.0, 340.0
    beta1 = 1.0 / (KB_KJ_PER_MOL_K * T1)
    beta2 = 1.0 / (KB_KJ_PER_MOL_K * T2)
    e1 = _synthetic_boltzmann_pe(rng, beta1, 200_000)
    e2 = _synthetic_boltzmann_pe(rng, beta2, 200_000)

    fit = shirts_slope_estimator_a(e1, T1, e2, T2)
    assert abs(fit.z) < 3.0, (
        f"Estimator A z-score too large on synthetic iid Boltzmann data: "
        f"slope={fit.slope:.5f} expected={fit.expected_slope:.5f} "
        f"se={fit.se:.5f} z={fit.z:.3f}"
    )


def test_shirts_estimator_b_recovers_known_slope_on_synthetic_boltzmann_samples():
    rng = np.random.default_rng(2025)
    T1, T2 = 280.0, 340.0
    beta1 = 1.0 / (KB_KJ_PER_MOL_K * T1)
    beta2 = 1.0 / (KB_KJ_PER_MOL_K * T2)
    e1 = _synthetic_boltzmann_pe(rng, beta1, 200_000)
    e2 = _synthetic_boltzmann_pe(rng, beta2, 200_000)

    fit = shirts_slope_estimator_b(e1, T1, e2, T2)
    assert abs(fit.z) < 3.0, (
        f"Estimator B z-score too large on synthetic iid Boltzmann data: "
        f"slope={fit.slope:.5f} expected={fit.expected_slope:.5f} "
        f"se={fit.se:.5f} z={fit.z:.3f}"
    )


def test_shirts_estimator_wrong_sign_is_rejected():
    """Sanity: feeding (e2,T2,e1,T1) SWAPPED must flip the expected sign but
    the (unswapped) FITTED slope stays put -- so a genuinely broken estimator
    (or a caller that mixes up the (e1,T1,e2,T2) order) shows up as a large
    |z|, not a silently-passing test. This guards the sign-discipline
    convention documented in the module docstring.
    """
    rng = np.random.default_rng(99)
    T1, T2 = 280.0, 340.0
    beta1 = 1.0 / (KB_KJ_PER_MOL_K * T1)
    beta2 = 1.0 / (KB_KJ_PER_MOL_K * T2)
    e1 = _synthetic_boltzmann_pe(rng, beta1, 100_000)
    e2 = _synthetic_boltzmann_pe(rng, beta2, 100_000)

    correct = shirts_slope_estimator_a(e1, T1, e2, T2)
    swapped_labels_only = shirts_slope_estimator_a(e1, T2, e2, T1)  # wrong T assignment
    assert abs(correct.z) < 3.0
    assert abs(swapped_labels_only.z) > 3.0, (
        "mislabelling T1/T2 should be caught by a large z-score, not pass silently"
    )


def test_thin_to_effective_keeps_most_of_an_iid_series():
    rng = np.random.default_rng(3)
    x = rng.normal(size=10_000)
    thinned = thin_to_effective(x)
    assert thinned.size > 0.5 * x.size  # g~1 for iid data -> stride 1


def test_parse_moves_csv_round_trip(tmp_path: Path):
    """A synthetic moves.csv round-trips through parse_moves_csv/world_series."""
    path = tmp_path / "synthetic.moves.csv"
    rows = [
        "round,replica,T,world,type,kick,accepted,PE,KE,fixman,H",
        "0,0,300,0,cartesian,0,1,-10.0,5.0,0.0,-5.0",
        "0,0,300,1,torsional,0,1,-10.0,2.0,-1.0,-9.0",
        "1,0,300,0,cartesian,0,1,-11.0,6.0,0.0,-5.0",
        "1,0,300,1,torsional,0,0,-11.0,3.0,-1.0,-9.0",
        "0,1,320,0,cartesian,0,1,-8.0,5.0,0.0,-3.0",
    ]
    path.write_text("\n".join(rows) + "\n")

    df = parse_moves_csv(path)
    assert list(df.columns) == MOVES_CSV_COLUMNS
    assert len(df) == 5

    pe_cart_r0 = world_series(df, "PE", world_type="cartesian", replica=0)
    np.testing.assert_allclose(pe_cart_r0, [-10.0, -11.0])

    pe_tors_r0 = world_series(df, "PE", world_type="torsional", replica=0)
    np.testing.assert_allclose(pe_tors_r0, [-10.0, -11.0])

    pe_r1 = world_series(df, "PE", replica=1)
    np.testing.assert_allclose(pe_r1, [-8.0])

    with pytest.raises(ValueError):
        world_series(df, "PE", world_type="does-not-exist")


def test_parse_moves_csv_missing_column_raises(tmp_path: Path):
    path = tmp_path / "bad.moves.csv"
    path.write_text("round,replica,T,world,type,kick,accepted,PE,KE,fixman\n0,0,300,0,cartesian,0,1,-1.0,0.0,0.0\n")
    with pytest.raises(ValueError):
        parse_moves_csv(path)
