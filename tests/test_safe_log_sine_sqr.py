import numpy as np
import robosample
from hypothesis import given
from hypothesis import strategies as st


# Ground truth
def reference_log_sin_sq(x):
    return 2.0 * np.log(np.abs(np.sin(x)))


# Invariants: no NaN / Inf
@given(
    st.floats(min_value=-1e-2, max_value=1e-2, allow_nan=False, allow_infinity=False)
)
def test_near_zero_stability(x):
    val = robosample.robo_bindings.safe_log_sine_sqr(x)

    assert np.isfinite(val)
    assert not np.isnan(val)
    assert not np.isinf(val)


# Agreement with ground truth away from 0
@given(st.floats(min_value=1e-2, max_value=10.0, allow_nan=False, allow_infinity=False))
def test_matches_reference(x):
    val = robosample.robo_bindings.safe_log_sine_sqr(x)
    ref = reference_log_sin_sq(x)

    assert np.isclose(val, ref, rtol=1e-6, atol=1e-8)


# Symmetry
@given(st.floats(min_value=1e-4, max_value=10.0))
def test_even_function(x):
    f1 = robosample.robo_bindings.safe_log_sine_sqr(x)
    f2 = robosample.robo_bindings.safe_log_sine_sqr(-x)

    assert np.isclose(f1, f2, rtol=1e-6, atol=1e-8)


def test_continuity_at_threshold():
    delta = 1e-6

    left = robosample.robo_bindings.safe_log_sine_sqr(delta * 0.999)
    right = robosample.robo_bindings.safe_log_sine_sqr(delta * 1.001)

    assert np.isclose(left, right, rtol=1e-3, atol=1e-6)


def test_matches_exact_outside():
    xs = np.logspace(-2, -1, 20)

    for x in xs:
        val = robosample.robo_bindings.safe_log_sine_sqr(x)
        ref = 2.0 * np.log(abs(np.sin(x)))

        assert np.isclose(val, ref, rtol=1e-6)


def test_local_smoothness():
    delta = 1e-6
    dx = 1e-8

    x1 = delta - dx
    x2 = delta + dx

    f1 = robosample.robo_bindings.safe_log_sine_sqr(x1)
    f2 = robosample.robo_bindings.safe_log_sine_sqr(x2)

    assert np.isclose(f1, f2, rtol=1e-3)


# Monotonicity near zero (as |x| increases, value increases)
def test_monotonicity_near_zero():
    xs = np.linspace(1e-6, 1e-2, 20)

    prev = robosample.robo_bindings.safe_log_sine_sqr(xs[0])
    for x in xs[1:]:
        curr = robosample.robo_bindings.safe_log_sine_sqr(x)
        assert curr > prev
        prev = curr


# Stress test with random sampling
@given(st.floats(min_value=-1e3, max_value=1e3, allow_nan=False, allow_infinity=False))
def test_wide_range_stress(x):
    val = robosample.robo_bindings.safe_log_sine_sqr(x)

    # Should never return NaN or inf
    assert np.isfinite(val)

    # Should not explode to absurd magnitudes
    assert abs(val) < 1e6
