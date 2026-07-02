"""Tier 1 -- configurational PE-distribution ladder (ala-dipeptide).

Spec: ``docs/specs/ensemble-validation/20-tier1-pe-distribution-ladder.md``
(foundations: ``docs/specs/ensemble-validation/00-foundations.md``).

Three rungs, all built on the SAME matched physics (FF/T/nonbonded/GBSA/no
constraints -- foundations Sec. 6 "matched physics"):

* **Rung A** -- native OpenMM Langevin (BAOAB/``LangevinMiddleIntegrator``)
  MD, the external anchor. Built via
  ``robosample.openmm_validation._build_reference_system`` (the same
  machinery ``test_openmm_potential_energy.py`` trusts) on the OpenMM
  ``Reference`` platform.
* **Rung B** -- Robosample ``add_cartesian_world`` alone.
* **Rung C** -- Robosample ``add_cartesian_world`` + ``add_robotic_world``
  (phi/psi torsional, Fixman ON), Gibbs-alternated.

Rung B and C's PE series are both the ``type=="cartesian"`` rows of
``moves.csv``: the Cartesian world's OpenMM MD relaxes ALL bond/angle/torsion
DOF every round in both rungs, so the ``cartesian``-world PE column always
reflects the SAME (full, flexible) DOF as Rung A -- the foundations Sec. 4
DOF-matching precondition for a histogram-EQUALITY test (framings #3 energy
leak and #4 native<->Robosample-OpenMM).

Cartesian-world acceptance mode (AlwaysAccept, not MetropolisHastings)
------------------------------------------------------------------------
``World::generateSample`` for a Cartesian world does NOT pull back the
device (OpenMM) kinetic energy after the MD trajectory (``World.cpp``:
"device kinetic energy not pulled back here", ``state_.energy.ke = 0.0``),
so its Metropolis test compares ``PE_new`` against ``PE_old`` alone, not the
full Hamiltonian ``H = PE+KE``. Empirically (verified interactively before
writing this suite) this makes ``MetropolisHastings`` on a Cartesian world
reject essentially every move (every fresh momentum draw's kinetic energy
that flows into potential energy during the trajectory reads as a "PE
increase" with no compensating KE decrease in the test), freezing the chain
at its start coordinates -- not a meaningful canonical sample. This is a
pre-existing implementation property of the Cartesian-world sampler, not a
theory question this spec touches, and the spec explicitly permits either
mode ("AlwaysAccept or MetropolisHastings"). Rung B/C therefore use
``AlwaysAccept``, i.e. trust the resample-then-integrate scheme (Andersen-
style thermostat + symplectic Verlet) to sample canonically without an
extra Metropolis correction -- the same principle Rung A's Langevin
integrator relies on (also no Metropolis step).

Sampling parameters (documented per CLAUDE.md test-gate requirement)
------------------------------------------------------------------------
Kept short but sufficient for the statistics below to have a discriminating
sample; longer runs only tighten the bands, they do not change what is being
tested. All three rungs sample every ``STRIDE``-equivalent unit of MD time
per logged point (Rung A: every 20 x 1fs steps between state queries; Rung
B/C cartesian world: every round == 20 x 1fs MD steps), so the between-sample
spacing is matched across rungs.
"""

from __future__ import annotations

import os
import pathlib

import numpy as np
import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))
_SLOW = os.environ.get("ROBOSAMPLE_SLOW_TESTS") is not None

# Same gating discipline as test_openmm_potential_energy.py: skip silently on
# a bare dev run when OpenMM/the compiled extension are missing, but HARD
# FAIL under the authoritative gate (ROBOSAMPLE_REQUIRE_OPENMM=1) so a
# missing dependency cannot silently shrink the test suite.
if _REQUIRE_OPENMM:
    import openmm as mm
    import openmm.unit as unit
    import robosample
else:
    pytest.importorskip("openmm")
    mm = pytest.importorskip("openmm")
    unit = pytest.importorskip("openmm.unit")
    robosample = pytest.importorskip("robosample")
from robosample import openmm_validation  # noqa: E402  (after importorskip)

import ensemble_stats as es  # noqa: E402  (sibling module in tests/)

# Every test in this module runs real MD/HMC sampling and SHALL NOT run in
# the default `nox -s tests` gate (CLAUDE.md "Slow-test gate").
pytestmark = pytest.mark.skipif(
    not _SLOW, reason="slow: sampling run; set ROBOSAMPLE_SLOW_TESTS=1"
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PRMTOP = REPO_ROOT / "examples" / "ala-dipeptide.prmtop"
RST7 = REPO_ROOT / "examples" / "ala-dipeptide.rst7"

# T_LOW doubles as BOTH the histogram-equality temperature (B==A, C==B) and
# the lower endpoint of the Shirts slope pair, so only 2 distinct
# temperatures (not 3) are needed in total (foundations Sec. 6 temperature-
# gap guidance: a ~20% gap gives a real slope signal while keeping the two
# histograms overlapping).
T_LOW = 300.0
T_HIGH = 360.0

# >=3 seeds so the Shirts slope check can be judged "consistently" across
# repeats (foundations Sec. 6 / Shirts checks.md line 134), not on a single
# noisy draw. Rung A/B/C all reuse the SAME seed value per repeat (fed to
# each simulator's own RNG) so "repeat" means an independent trajectory on
# every rung, not a shared one.
SEEDS = (11, 22, 33)

# ---- Rung A: native OpenMM Langevin MD (the external anchor) --------------
_NATIVE_BURN_IN_STEPS = 5000
_NATIVE_N_SAMPLES = 6000
_NATIVE_STRIDE_STEPS = 20
_NATIVE_DT_PS = 0.001  # 1 fs
_NATIVE_FRICTION_PER_PS = 1.0

# ---- Rung B: Robosample Cartesian-only world -------------------------------
_RUNG_B_EQUIL_ROUNDS = 200
_RUNG_B_PROD_ROUNDS = 6000
_RUNG_B_MDSTEPS = 20
_RUNG_B_DT_PS = 0.001  # 1 fs -- 20fs of MD per logged round, matching Rung A's stride

# ---- Rung C: Robosample Cartesian + robotic (torsional) composite --------
_RUNG_C_EQUIL_ROUNDS = 200
_RUNG_C_PROD_ROUNDS = 2000
_RUNG_C_CART_MDSTEPS = 20
_RUNG_C_CART_DT_PS = 0.001  # 1 fs, same as Rung B's cartesian world
_RUNG_C_TORS_MDSTEPS = 50
_RUNG_C_TORS_DT_PS = 0.002  # 2 fs


def _skip_if_missing_inputs():
    if _REQUIRE_OPENMM:
        assert PRMTOP.exists() and RST7.exists(), (
            f"example inputs not found: {PRMTOP} / {RST7} -- the authoritative "
            "gate (ROBOSAMPLE_REQUIRE_OPENMM=1) requires them to be present"
        )
    elif not PRMTOP.exists() or not RST7.exists():
        pytest.skip(f"example inputs not found: {PRMTOP} / {RST7}")


def _native_langevin_pe_series(topo_context, T: float, *, seed: int) -> np.ndarray:
    """Rung A: native OpenMM Langevin MD PE series at temperature T.

    Reuses ``openmm_validation._build_reference_system`` so the System
    (force field, GBSA, nonbonded method/cutoff) is IDENTICAL to what the
    Robosample rungs build from the same prmtop/rst7 -- foundations Sec. 6
    "matched physics".
    """
    inpcrd, system = openmm_validation._build_reference_system(
        topo_context, str(PRMTOP), str(RST7)
    )
    integrator = mm.LangevinMiddleIntegrator(
        T * unit.kelvin,
        _NATIVE_FRICTION_PER_PS / unit.picosecond,
        _NATIVE_DT_PS * unit.picoseconds,
    )
    integrator.setRandomNumberSeed(seed)
    platform = mm.Platform.getPlatformByName("Reference")
    mm_context = mm.Context(system, integrator, platform)
    mm_context.setPositions(inpcrd.positions)
    mm_context.setVelocitiesToTemperature(T * unit.kelvin, seed)

    mm_context.getIntegrator().step(_NATIVE_BURN_IN_STEPS)
    pes = np.empty(_NATIVE_N_SAMPLES)
    for i in range(_NATIVE_N_SAMPLES):
        mm_context.getIntegrator().step(_NATIVE_STRIDE_STEPS)
        state = mm_context.getState(getEnergy=True)
        pes[i] = state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)
    return pes


def _run_robosample_cartesian(base_name: str, T: float, *, seed: int) -> np.ndarray:
    """Rung B: Robosample Cartesian-only world PE series at temperature T."""
    context = robosample.Context(base_name, seed, robosample.AmberDihedralClassifier())
    context.load_amber(str(PRMTOP), str(RST7))
    context.set_enforce_periodic_box(False)
    context.add_cartesian_world().add_sampler(
        timeStep=_RUNG_B_DT_PS,
        mdSteps=_RUNG_B_MDSTEPS,
        acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
        use_nuts=False,
    )
    context.initialize([T])
    context.run_rex(_RUNG_B_EQUIL_ROUNDS, _RUNG_B_PROD_ROUNDS, 1, False)

    df = es.parse_moves_csv(base_name + ".moves.csv")
    return es.world_series(
        df, "PE", world_type="cartesian", replica=0, min_round=_RUNG_B_EQUIL_ROUNDS
    )


def _phi_psi_selection(context):
    df = context.standard_dihedral_bonds
    dihedrals = [
        robosample.DihedralType.PROTEIN_PHI.value,
        robosample.DihedralType.PROTEIN_PSI.value,
    ]
    bonds = df.loc[df["dihedral_type"].isin(dihedrals)]
    return context.build_flexibilities(bonds, robosample.rb.JointType.Torsion, False)


def _run_robosample_composite(base_name: str, T: float, *, seed: int) -> np.ndarray:
    """Rung C: Robosample Cartesian + robotic (phi/psi torsional, Fixman ON) composite.

    Returns the ``type=="cartesian"`` PE column -- the DOF-matched series
    (module docstring): the Cartesian world relaxes the full flexible system
    every round, in both this rung and Rung B.
    """
    context = robosample.Context(base_name, seed, robosample.AmberDihedralClassifier())
    context.load_amber(str(PRMTOP), str(RST7))
    context.set_enforce_periodic_box(False)
    # Cartesian world added FIRST so the Gibbs sweep runs it before the
    # torsional world every round (Context::runREX sweeps worlds_ in add
    # order) -- required for the "cartesian" PE column to reflect a fresh
    # full-DOF relaxation on every logged round.
    context.add_cartesian_world().add_sampler(
        timeStep=_RUNG_C_CART_DT_PS,
        mdSteps=_RUNG_C_CART_MDSTEPS,
        acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
        use_nuts=False,
    )
    sele = _phi_psi_selection(context)
    context.add_robotic_world(sele).add_sampler(
        timeStep=_RUNG_C_TORS_DT_PS,
        mdSteps=_RUNG_C_TORS_MDSTEPS,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        use_nuts=False,
        use_fixman=True,  # Fixman ON (foundations Sec. 3 binding convention)
    )
    context.initialize([T])
    context.run_rex(_RUNG_C_EQUIL_ROUNDS, _RUNG_C_PROD_ROUNDS, 1, False)

    df = es.parse_moves_csv(base_name + ".moves.csv")
    return es.world_series(
        df, "PE", world_type="cartesian", replica=0, min_round=_RUNG_C_EQUIL_ROUNDS
    )


@pytest.fixture(scope="module")
def pe_ladder_runs(tmp_path_factory):
    """runs[rung][seed][T] -> PE ndarray, for rung in {"A","B","C"}."""
    _skip_if_missing_inputs()
    workdir = tmp_path_factory.mktemp("pe_ladder")

    topo_context = robosample.Context(
        str(workdir / "topo_probe"), 0, robosample.AmberDihedralClassifier()
    )
    topo_context.load_amber(str(PRMTOP), str(RST7))
    topo_context.set_enforce_periodic_box(False)

    runs = {"A": {}, "B": {}, "C": {}}
    for seed in SEEDS:
        runs["A"][seed] = {
            T: _native_langevin_pe_series(topo_context, T, seed=seed)
            for T in (T_LOW, T_HIGH)
        }
        runs["B"][seed] = {
            T: _run_robosample_cartesian(
                str(workdir / f"rungB_s{seed}_T{int(T)}"), T, seed=seed
            )
            for T in (T_LOW, T_HIGH)
        }
        runs["C"][seed] = {
            T: _run_robosample_composite(
                str(workdir / f"rungC_s{seed}_T{int(T)}"), T, seed=seed
            )
            for T in (T_LOW, T_HIGH)
        }
    return runs


# ---------------------------------------------------------------------------
# B == A and C == B: DOF-matched PE-histogram equality (foundations Sec. 4)
# ---------------------------------------------------------------------------


def test_rung_b_matches_rung_a_pe_histogram(pe_ladder_runs):
    """INVARIANT (B==A): KS p>1e-3 at T_LOW, on the first seed's run."""
    a = pe_ladder_runs["A"][SEEDS[0]][T_LOW]
    b = pe_ladder_runs["B"][SEEDS[0]][T_LOW]
    stat, p, n1, n2 = es.ks_test_neff(a, b)
    assert p > 1e-3, (
        f"Rung B (Robosample Cartesian) PE histogram != Rung A (native OpenMM "
        f"Langevin) at T={T_LOW}: KS stat={stat:.4f} p={p:.3e} "
        f"(n_eff: A={n1}, B={n2})"
    )


def test_rung_c_matches_rung_b_pe_histogram(pe_ladder_runs):
    """INVARIANT (C==B): the energy-leak test (framing #3).

    A plausible-but-biased Fixman sign/normalization leaves B==A intact
    (Rung A/B never touch Fixman) but breaks C==B, because the composite's
    torsional-world moves would then bias the shared "cartesian" PE column
    away from the pure-Cartesian rung's distribution -- this is what
    discriminates the leak from a merely-shifted-but-still-canonical rung C
    (see also the per-rung slope check below on rung C).
    """
    b = pe_ladder_runs["B"][SEEDS[0]][T_LOW]
    c = pe_ladder_runs["C"][SEEDS[0]][T_LOW]
    stat, p, n1, n2 = es.ks_test_neff(b, c)
    assert p > 1e-3, (
        f"Rung C (Cartesian+GCHMC composite) PE histogram != Rung B (Cartesian "
        f"only) at T={T_LOW}: KS stat={stat:.4f} p={p:.3e} "
        f"(n_eff: B={n1}, C={n2})"
    )


# ---------------------------------------------------------------------------
# Moments: coarse pre-filter (foundations Sec. 6)
# ---------------------------------------------------------------------------


def test_moments_rung_b_matches_rung_a(pe_ladder_runs):
    a = pe_ladder_runs["A"][SEEDS[0]][T_LOW]
    b = pe_ladder_runs["B"][SEEDS[0]][T_LOW]
    ok, diff, tol = es.means_within_tolerance(a, b)
    assert ok, f"<PE> Rung B vs Rung A at T={T_LOW}: |diff|={diff:.4f} > 4*stderr={tol:.4f}"


def test_moments_rung_c_matches_rung_b(pe_ladder_runs):
    b = pe_ladder_runs["B"][SEEDS[0]][T_LOW]
    c = pe_ladder_runs["C"][SEEDS[0]][T_LOW]
    ok, diff, tol = es.means_within_tolerance(b, c)
    assert ok, f"<PE> Rung C vs Rung B at T={T_LOW}: |diff|={diff:.4f} > 4*stderr={tol:.4f}"


# ---------------------------------------------------------------------------
# Shirts PE-slope, every rung, consistently across seeds (foundations Sec. 5.1/6)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("rung", ["A", "B", "C"])
def test_shirts_slope_per_rung_consistent_across_seeds(pe_ladder_runs, rung):
    """LEMMA: per-rung Shirts slope == -(beta2-beta1) within 3*se(N_eff),
    consistently across seeds (mean |z| < 3 for both Estimator A and B).

    This validates canonical-ness of EACH rung independently of the B==A==C
    histogram-equality tests above (a rung could pass the slope test while
    disagreeing with its neighbours on the DOF-matched comparison, or vice
    versa -- the two are complementary, foundations Sec. 6 "choice of
    two-sample statistic").
    """
    z_a, z_b = [], []
    for seed in SEEDS:
        lo = pe_ladder_runs[rung][seed][T_LOW]
        hi = pe_ladder_runs[rung][seed][T_HIGH]
        fit_a = es.shirts_slope_estimator_a(lo, T_LOW, hi, T_HIGH)
        fit_b = es.shirts_slope_estimator_b(lo, T_LOW, hi, T_HIGH)
        z_a.append(fit_a.z)
        z_b.append(fit_b.z)

    mean_abs_z_a = float(np.mean(np.abs(z_a)))
    mean_abs_z_b = float(np.mean(np.abs(z_b)))
    assert mean_abs_z_a < 3.0, (
        f"rung {rung}: Estimator A z-scores across seeds {SEEDS} = "
        f"{[f'{z:.2f}' for z in z_a]} (mean|z|={mean_abs_z_a:.2f})"
    )
    assert mean_abs_z_b < 3.0, (
        f"rung {rung}: Estimator B z-scores across seeds {SEEDS} = "
        f"{[f'{z:.2f}' for z in z_b]} (mean|z|={mean_abs_z_b:.2f})"
    )
