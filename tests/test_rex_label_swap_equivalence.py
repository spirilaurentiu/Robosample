"""INVARIANT-EQUIV -- label-swap REMC vs legacy coordinate-swap REMC.

Spec: ``docs/specs/replica-exchange-nonequilibrium-work.md``, "Port-target
baseline" ("the two are provably equivalent for pure-temperature REMC (same
permutation distribution over states)") and the Stage 1 coder checkpoint
(2-replica smoke found the equivalence is not merely statistical but EXACT,
once ``RunREX`` visits the shared worlds in the same per-round order as the
legacy driver -- see ``Context::RunREX``'s Gibbs-sweep comment).

Regression scope (>=3-replica, reviewer Should-fix follow-up on Stage 1): the
2-replica smoke run the coder used interactively can only ever exercise ONE
``prepareExchangePairs`` parity class -- with T==2 there is exactly one
possible pair, (0,1), regardless of the round's oddity. It therefore cannot
catch the B7 parity bug the spec calls out explicitly (revision 2, reviewer
Should-fix R5): if the exchange-round parity were ever taken from the wrong
counter, one whole pairing class (e.g. every ``(1,2),(3,4),...`` round) could
be silently skipped without this failing. This test uses T==5 thermodynamic
states so BOTH parity classes are structurally distinct and non-empty:

  * even parity (``exchangeRound_ % 2 == 0``): pairs (0,1), (2,3)
  * odd  parity (``exchangeRound_ % 2 == 1``): pairs (1,2), (3,4)

and asserts swaps were attempted in pairs from BOTH classes before trusting
the equivalence result.

Run: pytest -q tests/test_rex_label_swap_equivalence.py
Slow-test gate: set ROBOSAMPLE_SLOW_TESTS=1 (real MD/HMC sampling; skipped by
default under ``nox -s tests``, matching test_ensemble_pe_ladder.py).
"""

from __future__ import annotations

import csv
import os
import pathlib

import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))
_SLOW = os.environ.get("ROBOSAMPLE_SLOW_TESTS") is not None

if _REQUIRE_OPENMM:
    import robosample
else:
    robosample = pytest.importorskip("robosample")

pytestmark = pytest.mark.skipif(
    not _SLOW, reason="slow: sampling run; set ROBOSAMPLE_SLOW_TESTS=1"
)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
PRMTOP = REPO_ROOT / "examples" / "ala-dipeptide.prmtop"
RST7 = REPO_ROOT / "examples" / "ala-dipeptide.rst7"

# T=5 thermodynamic states (>=3 required; 5 gives 2 pairs in EACH parity
# class -- see module docstring). Modest gaps keep acceptance well away from
# both 0 and 1 on this tiny (32-atom) system.
TEMPERATURES = [300.0, 310.0, 320.0, 330.0, 340.0]
SEED = 424242
EQUIL_ROUNDS = 5
PROD_ROUNDS = 40
WRITE_FREQ = 1

# Both classes need >=2 EXECUTED mixes to fire (swapEvery=1 default => one
# mix per round, alternating parity every round) -- 45 total rounds is ample.


def _skip_if_missing_inputs():
    if _REQUIRE_OPENMM:
        assert PRMTOP.exists() and RST7.exists(), (
            f"example inputs not found: {PRMTOP} / {RST7} -- the authoritative "
            "gate (ROBOSAMPLE_REQUIRE_OPENMM=1) requires them to be present"
        )
    elif not PRMTOP.exists() or not RST7.exists():
        pytest.skip(f"example inputs not found: {PRMTOP} / {RST7}")


def _build_context(base_name: str, seed: int) -> "robosample.Context":
    context = robosample.Context(base_name, seed, robosample.AmberDihedralClassifier())
    context.load_amber(str(PRMTOP), str(RST7))
    context.set_enforce_periodic_box(False)

    dihedrals = [
        robosample.DihedralType.PROTEIN_PHI.value,
        robosample.DihedralType.PROTEIN_PSI.value,
    ]
    bonds = context.standard_dihedral_bonds.loc[
        context.standard_dihedral_bonds["dihedral_type"].isin(dihedrals)
    ]
    sele = context.build_flexibilities(bonds, robosample.rb.JointType.Torsion, False)
    context.add_robotic_world(sele).add_sampler(
        timeStep=0.002,
        mdSteps=20,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        use_nuts=False,
        use_fixman=True,
    )
    context.initialize(TEMPERATURES)
    return context


def _read_state_pe(base_name: str, state_idx: int) -> list[float]:
    path = f"{base_name}.{state_idx}.csv"
    pes = []
    with open(path, newline="") as f:
        for row in csv.reader(f):
            if len(row) != 4:
                continue
            pes.append(float(row[3]))
    return pes


@pytest.fixture(scope="module")
def rex_runs(tmp_path_factory):
    _skip_if_missing_inputs()
    workdir = tmp_path_factory.mktemp("rex_equiv")

    legacy_base = str(workdir / "legacy")
    ctx_legacy = _build_context(legacy_base, SEED)
    ctx_legacy.run_rex(EQUIL_ROUNDS, PROD_ROUNDS, WRITE_FREQ, False)

    label_base = str(workdir / "label")
    ctx_label = _build_context(label_base, SEED)
    ctx_label.set_replica_mixing_scheme(robosample.rb.ReplicaMixingScheme.Neighboring)
    ctx_label.set_swap_every(1)
    ctx_label.run_rex_label_swap(robosample.rb.RunType.REMC, EQUIL_ROUNDS, PROD_ROUNDS, WRITE_FREQ, False)

    attempted = ctx_label.attempted_swaps_matrix()
    return {
        "legacy_base": legacy_base,
        "label_base": label_base,
        "attempted": attempted,
    }


def test_both_b7_parity_classes_were_attempted(rex_runs):
    """Guards the regression this test exists for (B7 revision-2, R5): a
    dedicated exchange-round counter, not ``mixi % 2``/round parity aliasing,
    SHALL keep BOTH neighbour-pairing classes reachable. If parity ever
    froze on one class, (1,2)/(3,4) (or (0,1)/(2,3)) would read zero here.
    """
    attempted = rex_runs["attempted"]
    even_class_pairs = [(0, 1), (2, 3)]
    odd_class_pairs = [(1, 2), (3, 4)]
    for i, j in even_class_pairs + odd_class_pairs:
        assert attempted[i][j] > 0, (
            f"pair ({i},{j}) was never attempted -- one B7 parity class is "
            f"unreachable (full matrix: {attempted})"
        )


def test_per_state_pe_series_match_legacy_oracle(rex_runs):
    """INVARIANT-EQUIV: label-swap RunREX and legacy coordinate-swap runREX
    SHALL produce the same per-thermodynamic-state PE series. This checks
    exact (not merely statistical) equality -- the Stage 1 checkpoint's
    2-replica smoke found the two mechanisms are the SAME Markov chain
    realization bit-for-bit once RunREX iterates the Gibbs sweep by
    thermodynamic-state index (Context::RunREX doc comment); this test
    extends that check to >=3 replicas and both B7 parity classes.
    """
    T = len(TEMPERATURES)
    for k in range(T):
        pe_legacy = _read_state_pe(rex_runs["legacy_base"], k)
        pe_label = _read_state_pe(rex_runs["label_base"], k)
        assert len(pe_legacy) == len(pe_label) == PROD_ROUNDS, (
            f"state {k}: unexpected sample count (legacy={len(pe_legacy)}, "
            f"label={len(pe_label)}, expected {PROD_ROUNDS})"
        )
        for round_idx, (a, b) in enumerate(zip(pe_legacy, pe_label)):
            assert a == pytest.approx(b, abs=1e-6), (
                f"state {k} round {round_idx}: legacy PE={a!r} != "
                f"label-swap PE={b!r}"
            )
