"""V8 - exact acceptance algebra for the RENE/RENEMC replica-exchange swap.

Reviewer reproducer for docs/specs/replica-exchange-nonequilibrium-work.md
(Derivation sketch, B6, D2, V8). Noise-free spec oracle: it re-implements the
spec's log-alpha formula on a one-body analytic system and checks it against
first-principles detailed balance for the deterministic involution swap.

One-body system: each replica is a single bond, a vector x in R^3 with
U(x) = 0.5 k (|x| - r0)^2. The BAT scaling drives the bond isotropically,
x -> s x, whose exact Cartesian log-Jacobian is 3 ln s (D6: a single isotropic
bond scale gives |dx'/dx| = s^3). Closed form, no MD, no RNG.

What this file establishes (all three verified by re-derivation, see the
reviewer findings that ship with it):

  A. The spec's log-alpha WITH the correct beta-assignment (x^0 at source beta,
     x^tau at target beta) satisfies the true detailed-balance equation
     pi(z) a(z->Tz) = pi(Tz) a(Tz->z); the swapped beta-assignment does NOT.
     -> guards the beta-assignment.

  B. On the symmetric s / (1/s) two-replica involution the JOINT Jacobian is
     identically 1 (s^3 * s^-3), so a flipped Jacobian sign or a dropped ln s is
     INVISIBLE to the swap detailed-balance check. This contradicts spec V8's
     claim (spec line 589) that "a flipped Jacobian sign or a missing ln s SHALL
     make this fail". V8 as worded cannot fail on the Jacobian bug it names.
     -> documents the FINDING: V8 is not a Jacobian oracle.

  C. The Jacobian sign/magnitude is only observable where |det dT/dz| != 1, e.g.
     a single-replica NCMC radial drive (Nilmeier eq 28,
     A = min{1, e^{-beta[U(x_tau)-U(x_0)]} (r_tau/r_0)^... }, here full 3D
     isotropic so the factor is s^3). There a flipped sign or a dropped ln s
     diverges from the analytic acceptance.
     -> the guard the spec must actually use for the Jacobian (its V5).

Run: pytest -q tests/test_rex_swap_acceptance_algebra.py
"""

import numpy as np
import pytest

K_BOND = 3.0
R0 = 1.0


def U(x):
    r = np.linalg.norm(x)
    return 0.5 * K_BOND * (r - R0) ** 2


def lnjac_isotropic(s):
    """Exact Cartesian log-Jacobian of x -> s x in R^3 (D6: 3 ln s)."""
    return 3.0 * np.log(s)


def swap_logalpha(cold_cfg, hot_cfg, beta_C, beta_H, s,
                  jac_sign=+1.0, use_lns=True, swap_betas=False):
    """Spec WTerm = -(Work_X + Work_Y) for the paired drive + label swap.

    X = replica in the cold state, driven by s toward hot.
    Y = replica in the hot state, driven by 1/s toward cold.
    mean = 0 so M_s(q) = s q is an exact involution with M_{1/s}.

    jac_sign/use_lns/swap_betas inject the bugs the oracle must catch.
    """
    x_Xtau = s * cold_cfg          # post-scale config of X (x' == x^tau, no MD)
    x_Ytau = (1.0 / s) * hot_cfg   # post-scale config of Y
    lnJac_X = (lnjac_isotropic(s) if use_lns else 0.0) * jac_sign
    lnJac_Y = (lnjac_isotropic(1.0 / s) if use_lns else 0.0) * jac_sign

    if not swap_betas:
        # correct: x^tau at TARGET beta, x^0 at SOURCE beta
        Work_X = beta_H * U(x_Xtau) - beta_C * U(cold_cfg) - lnJac_X
        Work_Y = beta_C * U(x_Ytau) - beta_H * U(hot_cfg) - lnJac_Y
    else:
        # BUG: x^tau at source beta, x^0 at target beta
        Work_X = beta_C * U(x_Xtau) - beta_H * U(cold_cfg) - lnJac_X
        Work_Y = beta_H * U(x_Ytau) - beta_C * U(hot_cfg) - lnJac_Y

    WTerm = -(Work_X + Work_Y)
    return WTerm, (x_Xtau, x_Ytau)


def log_pi(cold_cfg, hot_cfg, beta_C, beta_H):
    """Joint target: cold state holds cold_cfg, hot state holds hot_cfg."""
    return -beta_C * U(cold_cfg) - beta_H * U(hot_cfg)


def _accept(logalpha):
    return min(1.0, np.exp(logalpha))


# --------------------------------------------------------------------------- #
# A. beta-assignment is guarded by the true detailed-balance equation.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("s", [1.15, 1.4, 0.8])
@pytest.mark.parametrize("betas", [(1.0, 0.5), (2.0, 1.3)])
def test_beta_assignment_satisfies_detailed_balance(s, betas):
    beta_C, beta_H = betas
    rng = np.random.default_rng(1)
    coldc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])
    hotc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])

    la, (xXt, xYt) = swap_logalpha(coldc, hotc, beta_C, beta_H, s)
    # Tz: cold state now holds the driven Y config, hot state holds driven X.
    la_rev, _ = swap_logalpha(xYt, xXt, beta_C, beta_H, s)

    lhs = np.exp(log_pi(coldc, hotc, beta_C, beta_H)) * _accept(la)
    rhs = np.exp(log_pi(xYt, xXt, beta_C, beta_H)) * _accept(la_rev)
    assert np.isclose(lhs, rhs, rtol=1e-12, atol=1e-14)


@pytest.mark.parametrize("s", [1.15, 1.4])
def test_swapped_beta_assignment_breaks_detailed_balance(s):
    """The x^0<->x^tau beta swap must be caught by the DB equation."""
    beta_C, beta_H = 1.0, 0.5
    rng = np.random.default_rng(2)
    coldc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])
    hotc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])

    la, (xXt, xYt) = swap_logalpha(coldc, hotc, beta_C, beta_H, s, swap_betas=True)
    la_rev, _ = swap_logalpha(xYt, xXt, beta_C, beta_H, s, swap_betas=True)

    lhs = np.exp(log_pi(coldc, hotc, beta_C, beta_H)) * _accept(la)
    rhs = np.exp(log_pi(xYt, xXt, beta_C, beta_H)) * _accept(la_rev)
    assert not np.isclose(lhs, rhs, rtol=1e-9, atol=1e-12)


# --------------------------------------------------------------------------- #
# B. FINDING: the symmetric swap CANNOT see the Jacobian (joint |det| == 1).
#    These asserts pin the gap so nobody trusts V8 to guard D6/F7.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("s", [1.15, 1.4])
def test_symmetric_swap_is_blind_to_jacobian_sign_flip(s):
    beta_C, beta_H = 1.0, 0.5
    rng = np.random.default_rng(3)
    coldc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])
    hotc = rng.normal(size=3) + np.array([1.0, 0.0, 0.0])

    def db_gap(**kw):
        la, (xXt, xYt) = swap_logalpha(coldc, hotc, beta_C, beta_H, s, **kw)
        la_rev, _ = swap_logalpha(xYt, xXt, beta_C, beta_H, s, **kw)
        lhs = np.exp(log_pi(coldc, hotc, beta_C, beta_H)) * _accept(la)
        rhs = np.exp(log_pi(xYt, xXt, beta_C, beta_H)) * _accept(la_rev)
        return abs(lhs - rhs)

    # Correct, sign-flipped, and ln-s-dropped all satisfy DB equally: the swap
    # move is blind to the Jacobian because s^3 * s^-3 == 1.
    assert db_gap() < 1e-12
    assert db_gap(jac_sign=-1.0) < 1e-12
    assert db_gap(use_lns=False) < 1e-12


# --------------------------------------------------------------------------- #
# C. The Jacobian is observable only where |det dT/dz| != 1: single-body NCMC.
#    This is the oracle the Jacobian bug (D6/F7) actually needs.
# --------------------------------------------------------------------------- #

@pytest.mark.parametrize("s", [1.15, 1.4, 0.8])
def test_single_body_ncmc_acceptance_needs_correct_jacobian(s):
    beta = 1.0
    rng = np.random.default_rng(4)
    x0 = rng.normal(size=3) + np.array([1.5, 0.0, 0.0])
    xtau = s * x0

    # Nilmeier eq 28 for a deterministic radial (here 3D isotropic) drive.
    analytic = -beta * (U(xtau) - U(x0)) + lnjac_isotropic(s)

    correct = -beta * (U(xtau) - U(x0)) + lnjac_isotropic(s)
    flipped = -beta * (U(xtau) - U(x0)) - lnjac_isotropic(s)
    dropped = -beta * (U(xtau) - U(x0))

    assert np.isclose(correct, analytic, rtol=1e-12, atol=1e-14)
    assert abs(flipped - analytic) > 1e-6   # flipped sign is caught here
    assert abs(dropped - analytic) > 1e-6   # missing ln s is caught here
