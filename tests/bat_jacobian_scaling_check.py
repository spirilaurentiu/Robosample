"""V5 analytic Jacobian oracle for the RENE BAT-scaling map (spec D6 / F7).

Noise-free discriminator for the corrected forward Cartesian log-Jacobian of the
deterministic scaling map  x0 -->(Cart->BAT) q0 -->(scale) q_tau -->(BAT->Cart) x'.

Spec claim (docs/specs/replica-exchange-nonequilibrium-work.md, D6):

    lnJac = ( J(x') - J(x0) ) + N_scaled * ln s ,   J(x) = sum_bodies [ 2 ln r + ln sin(theta) ]

Original (HMCSampler.cpp:659):  bendStretchJacobianDetLog = J_ini + J_scale - J_fin
  = ( J(x0) - J(x') ) + J_scale   -- inverted (J_fin - J_ini) sign + bogus J_scale.

This file re-derives the volume element and checks the corrected formula against a
finite-difference of the actual Cartesian map. It has no engine dependency; it is the
reference the ported getDistortJacobianDetLog() must reproduce (V5).

Run:  python3 tests/bat_jacobian_scaling_check.py     (also usable under pytest)
"""

import math
import numpy as np


# ---- spherical BAT <-> Cartesian for one atom about a fixed parent frame -------
# A = P + r * (sin th cos ph, sin th sin ph, cos th).  |dA/d(r,th,ph)| = r^2 sin th.
# So the per-body volume element J = ln(r^2 sin th) = 2 ln r + ln sin th (spec, CONFIRMED).

P = np.array([0.30, -0.20, 0.15])  # fixed parent origin (arbitrary, nonzero)


def bat_to_cart(r, th, ph):
    return P + r * np.array([math.sin(th) * math.cos(ph),
                             math.sin(th) * math.sin(ph),
                             math.cos(th)])


def cart_to_bat(x):
    v = x - P
    r = np.linalg.norm(v)
    th = math.acos(v[2] / r)
    ph = math.atan2(v[1], v[0])
    return r, th, ph


def J_of(x):
    """Geometric BAT volume-element log:  2 ln r + ln sin theta."""
    r, th, _ = cart_to_bat(x)
    return 2.0 * math.log(r) + math.log(math.sin(th))


def fd_log_jac(scaling_map, x0, h=1e-6):
    """log|det d x' / d x0| by central finite difference of the Cartesian map."""
    n = len(x0)
    Jm = np.zeros((n, n))
    for k in range(n):
        e = np.zeros(n); e[k] = h
        Jm[:, k] = (scaling_map(x0 + e) - scaling_map(x0 - e)) / (2.0 * h)
    sign, logdet = np.linalg.slogdet(Jm)
    assert sign > 0, "map orientation reversed"
    return logdet


# ---- the two elementary scaling maps -------------------------------------------
def bond_scale_map(s):
    def f(x):
        r, th, ph = cart_to_bat(x)
        return bat_to_cart(s * r, th, ph)      # scale bond only
    return f


def angle_scale_map(s):
    def f(x):
        r, th, ph = cart_to_bat(x)
        return bat_to_cart(r, s * th, ph)      # scale angle only (mean 0)
    return f


# ---- corrected vs original closed forms ----------------------------------------
def corrected_lnJac(x0, xt, n_scaled, s):
    return (J_of(xt) - J_of(x0)) + n_scaled * math.log(s)


def original_lnJac_bond(x0, xt, s, drop_double_count):
    """Reproduce HMCSampler.cpp J_ini + J_scale - J_fin for a single bond scale.
    J_ini - J_fin = -(J(x')-J(x0)) = -2 ln s.
    J_scale for a bond is accumulated at :590 (guarded) AND :624 (unconditional) -> 2 ln s;
    with the F8 double-count removed it is ln s."""
    j_ini_minus_fin = J_of(x0) - J_of(xt)          # = -2 ln s
    j_scale = (1.0 if drop_double_count else 2.0) * math.log(s)
    return j_ini_minus_fin + j_scale


def _check(name, got, want, tol=1e-6):
    ok = abs(got - want) < tol
    print(f"  {name:34s} got={got:+.8f} want={want:+.8f} {'OK' if ok else 'FAIL'}")
    return ok


def run():
    ok = True
    r0, th0, ph0 = 1.30, 0.9, 0.6
    x0 = bat_to_cart(r0, th0, ph0)

    for s in (1.25, 0.80, 1.5, 1.0):
        print(f"\n[s = {s}]")

        # --- single isotropic BOND scale: true Cartesian |dx'/dx0| = s^3 -> 3 ln s ---
        f = bond_scale_map(s)
        xt = f(x0)
        fd = fd_log_jac(f, x0)
        ok &= _check("FD Cartesian log-jac (bond)", fd, 3.0 * math.log(s))
        ok &= _check("corrected formula (bond)",
                     corrected_lnJac(x0, xt, n_scaled=1, s=s), 3.0 * math.log(s))
        # original is wrong: 0 (double-counted) or -ln s (F8-fixed); neither is 3 ln s
        orig0 = original_lnJac_bond(x0, xt, s, drop_double_count=False)
        origF8 = original_lnJac_bond(x0, xt, s, drop_double_count=True)
        ok &= _check("original J_ini+J_scale-J_fin", orig0, 0.0)
        ok &= _check("original after F8 fix", origF8, -math.log(s))
        if s != 1.0:
            assert abs(orig0 - 3 * math.log(s)) > 1e-3, "original must differ from truth"
            assert abs(origF8 - 3 * math.log(s)) > 1e-3, "F8-fixed original still biased"

        # --- single ANGLE scale: ln(sin th'/sin th) + ln s ---
        g = angle_scale_map(s)
        xt = g(x0)
        thp = s * th0
        want_angle = math.log(math.sin(thp) / math.sin(th0)) + math.log(s)
        ok &= _check("FD Cartesian log-jac (angle)", fd_log_jac(g, x0), want_angle)
        ok &= _check("corrected formula (angle)",
                     corrected_lnJac(x0, xt, n_scaled=1, s=s), want_angle)

    print("\nRESULT:", "PASS" if ok else "FAIL")
    return ok


def test_bat_jacobian_scaling():
    assert run()


if __name__ == "__main__":
    import sys
    sys.exit(0 if run() else 1)
