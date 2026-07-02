"""Shared helpers for the Fixman idealized-chain / proline-realism pytest tiers
(docs/specs/fixman-idealized-chains-validation.md Tier 1 / Tier 2).

Python full-stack surface (spec Sec. 3.3): drives the shipped
``examples/fixman/*`` prmtops through the ``Context.load_amber`` -> world ->
``run_rex`` path (precedent ``python/robosample/run.py``), histogramming
torsions read back from the DCD trajectory the sampler writes (mdtraj) -- the
ONLY per-round position readback the Python API exposes. ``SystemTopology.
atoms_x`` is populated once at ``load_amber`` time and never updated per round
(``Context::runREX`` keeps positions in a private ``replicaCoords_``; see
``docs/specs/ensemble-validation/30-tier2-torsion-conformational.md``
"Torsion extraction" for the same DCD+mdtraj convention this module follows).
"""

from __future__ import annotations

import math
import os

import numpy as np
import scipy.stats

import robosample


def slow_tests_enabled() -> bool:
    return bool(os.environ.get("ROBOSAMPLE_SLOW_TESTS"))


def rotatable_bonds(context) -> list[tuple[int, int]]:
    """Every non-ring-closing bond whose BOTH endpoints carry >= 2 bonds --
    ``Context::buildFlexibilities``'s own C++ "rotatable" rule (Context.cpp
    ~:143-144), replicated here in Python.

    NOTE (pre-existing, unrelated bug, flagged not fixed -- out of scope for
    this spec): ``Context.build_flexibilities(bonds=None, ...)`` is documented
    (both the Python and C++ docstrings) to mean "every eligible bond", but
    ``context.py``'s wrapper coerces ``None`` to ``[]`` before the C++ call,
    which pybind11 turns into a POPULATED-but-EMPTY ``std::optional``, so
    ``Context::buildFlexibilities`` sees ``bonds.has_value()==true`` with an
    empty ``want`` set and flexes NOTHING (confirmed empirically:
    ``world.n_dof == 0`` after ``build_flexibilities(None, ...)`` on
    ``examples/fixman/c4``). Passing an EXPLICIT bond list (this function's
    result) sidesteps the bug without touching ``context.py``.
    """
    st = context.system_topology
    bonds_i, bonds_j = list(st.bonds_i), list(st.bonds_j)
    ring = list(st.bonds_ring_closing)
    deg = list(st.atoms_num_bonds_involved)
    return [
        (i, j)
        for i, j, r in zip(bonds_i, bonds_j, ring)
        if (not r) and deg[i] >= 2 and deg[j] >= 2
    ]


def build_torsional_context(
    base_name: str,
    prmtop: str,
    rst7: str,
    seed: int,
    use_fixman: bool,
    bonds: list[tuple[int, int]] | None = None,
    timestep: float = 0.002,
    md_steps: int = 100,
):
    """A torsional (internal-coordinate) world: TORSIONAL (use_fixman=False) or
    FIXMAN (use_fixman=True) per spec Sec. 2. Returns (context, world, bonds_used).

    Mixing default (timestep=0.002 ps, md_steps=100): a torsional world freezes
    bond/angle vibrations (rigid bodies), so a 2 fs step is stable even under a
    real force field, and 100 MD steps/round drops the torsion IACT from ~360
    (at the old 20/0.001) to ~10 -- WITHOUT this the per-round samples are ~30x
    autocorrelated and no goodness-of-fit test at a fixed alpha is valid (see
    ``thin_to_independent``).
    """
    ctx = robosample.Context(base_name, seed, robosample.AmberDihedralClassifier())
    ctx.load_amber(str(prmtop), str(rst7), use_gbsa_obc2=False)
    sel_bonds = rotatable_bonds(ctx) if bonds is None else bonds
    sele = ctx.build_flexibilities(sel_bonds, robosample.rb.JointType.Torsion, False)
    world = ctx.add_robotic_world(sele)
    world.add_sampler(
        timeStep=timestep,
        mdSteps=md_steps,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        use_nuts=False,
        use_fixman=use_fixman,
    )
    return ctx, world, sel_bonds


def build_cartesian_context(
    base_name: str,
    prmtop: str,
    rst7: str,
    seed: int,
    timestep: float = 0.0005,
    md_steps: int = 50,
):
    """A Cartesian (OpenMM-MD) world: the FLEXIBLE reference ensemble."""
    ctx = robosample.Context(base_name, seed, robosample.AmberDihedralClassifier())
    ctx.load_amber(str(prmtop), str(rst7), use_gbsa_obc2=False)
    world = ctx.add_cartesian_world()
    world.add_sampler(
        timeStep=timestep,
        mdSteps=md_steps,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        use_nuts=False,
    )
    return ctx


def run_and_load_dihedrals(
    ctx,
    base_name: str,
    prmtop: str,
    equil_rounds: int,
    prod_rounds: int,
    write_freq: int,
    atom_quads: list[list[int]],
) -> np.ndarray:
    """Run equil+prod rounds (single replica, T=300K) and return an
    (n_frames, len(atom_quads)) array of dihedral angles [rad], read back from
    the DCD trajectory ``{base_name}.0.dcd`` the sampler writes.
    """
    import mdtraj

    ctx.initialize([300])
    ctx.run_rex(equil_rounds, prod_rounds, write_freq, False)
    traj = mdtraj.load(f"{base_name}.0.dcd", top=str(prmtop))
    return mdtraj.compute_dihedrals(traj, atom_quads)


# ---------------------------------------------------------------------------
#  Statistics: chi-square goodness-of-fit (Tier 1) and Hellinger distance
#  (Tier 2), matching StatTest.hpp's conventions (alpha = 1e-4).
# ---------------------------------------------------------------------------


def chi_square_critical(dof: int, alpha: float = 1e-4) -> float:
    return float(scipy.stats.chi2.isf(alpha, dof))


def chi_square_statistic(observed: np.ndarray, expected: np.ndarray) -> float:
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    mask = expected > 0
    return float(np.sum((observed[mask] - expected[mask]) ** 2 / expected[mask]))


def histogram_counts(
    values: np.ndarray, nbins: int, lo: float = -math.pi, hi: float = math.pi
):
    """Wrap values into [lo, hi) then histogram -- returns (counts, edges)."""
    span = hi - lo
    wrapped = lo + np.mod(np.asarray(values, dtype=float) - lo, span)
    counts, edges = np.histogram(wrapped, bins=nbins, range=(lo, hi))
    return counts, edges


def expected_counts_from_density(density_fn, edges: np.ndarray, total_n: int, quad_pts: int = 64):
    """Expected COUNTS per bin for an unnormalized density (evaluated by
    midpoint quadrature within each bin, then rescaled so the total equals
    ``total_n``) -- the same "expectedFromWeights" idea as
    ``StatTest.hpp::expectedFromWeights``, generalized to sub-bin quadrature so
    a sharply-peaked density (e.g. eq:24 near its extrema) is still integrated
    accurately at a modest bin count.
    """
    weights = np.empty(len(edges) - 1, dtype=float)
    for b in range(len(edges) - 1):
        lo, hi = edges[b], edges[b + 1]
        xs = lo + (np.arange(quad_pts) + 0.5) * (hi - lo) / quad_pts
        weights[b] = float(np.mean([density_fn(x) for x in xs]))
    total_w = weights.sum()
    return weights * (total_n / total_w)


def hellinger_distance(p: np.ndarray, q: np.ndarray) -> float:
    """Hellinger distance between two (not necessarily normalized) histograms
    over the SAME binning: H(p,q) = sqrt(1/2 sum (sqrt(p_i)-sqrt(q_i))^2), p,q
    normalized to sum 1 first.
    """
    p = np.asarray(p, dtype=float)
    q = np.asarray(q, dtype=float)
    p = p / p.sum()
    q = q / q.sum()
    return float(np.sqrt(0.5 * np.sum((np.sqrt(p) - np.sqrt(q)) ** 2)))


# ---------------------------------------------------------------------------
#  Autocorrelation handling (REQUIRED before any goodness-of-fit / distance on
#  HMC output). A chi-square GOF test at a fixed alpha, and Hellinger thresholds,
#  are only valid on APPROXIMATELY INDEPENDENT samples. Robosample's MD-HMC
#  produces a torsion series with a long integrated autocorrelation time (IACT):
#  at mdSteps=20/ts=0.001 the C4 torsion has tau ~ 360, so ~1e4 per-round frames
#  carry only ~30 independent draws. Treating them as independent inflates chi2
#  into the thousands and rejects EVERY hypothesis (including the true one) -- the
#  original Tier-1/Tier-2 failure mode. These helpers estimate the IACT and thin
#  the series to independent draws so the downstream statistic is meaningful.
# ---------------------------------------------------------------------------


def integrated_autocorr_time(theta: np.ndarray, max_lag: int = 2000) -> float:
    """IACT of a PERIODIC series via its unit-vector (cos,sin) embedding (so the
    +-pi wrap does not create spurious jumps), summed with Sokal's automatic
    window (stop when the window M >= 6*tau). Returns >= 1.0."""
    x = np.column_stack([np.cos(theta), np.sin(theta)]).astype(float)
    x -= x.mean(axis=0)
    n = len(x)
    var = float(np.sum(x * x) / n)
    if n < 4 or var <= 0.0:
        return 1.0
    acf = [1.0]
    for lag in range(1, min(max_lag, n // 2)):
        acf.append(float(np.sum(x[: n - lag] * x[lag:]) / n / var))
    acf = np.asarray(acf)
    tau = 1.0
    for m in range(1, len(acf)):
        tau = 1.0 + 2.0 * float(np.sum(acf[1 : m + 1]))
        if m >= 6.0 * tau:
            break
    return max(tau, 1.0)


def thin_to_independent(values: np.ndarray):
    """Subsample a 1D torsion series by ceil(IACT) so residual autocorrelation is
    negligible. Returns (thinned_values, tau, stride). Histogram GOF / Hellinger
    thresholds should be applied to the thinned series, and effective N is
    len(thinned)."""
    values = np.asarray(values, dtype=float)
    tau = integrated_autocorr_time(values)
    stride = max(1, int(math.ceil(tau)))
    return values[::stride], tau, stride


def flat_noise_floor(n_eff: int, nbins: int) -> float:
    """Expected Hellinger-to-flat of a TRULY flat distribution sampled with n_eff
    independent draws into nbins bins (finite-sample noise): ~ sqrt((nbins-1)/(8
    n_eff)). A measured H below ~2x this is 'flat within noise'; well above it is
    a real deviation. Used to size flatness assertions without a magic constant."""
    n_eff = max(int(n_eff), 1)
    return float(math.sqrt((nbins - 1) / (8.0 * n_eff)))
