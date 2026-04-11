import math
import struct

import mpmath as mp
import numpy as np
import pytest
import robosample
from hypothesis import given, settings
from hypothesis import strategies as st

# ============================================================
# Helpers
# ============================================================


def reference_lse2(a, b, prec=100):
    mp.mp.dps = prec

    if math.isinf(a) or math.isinf(b):
        if a == math.inf or b == math.inf:
            return math.inf
        if a == -math.inf and b == -math.inf:
            return -math.inf

    return float(mp.log(mp.e**a + mp.e**b))


def float_to_int(x):
    """Bit-level representation for ULP distance."""
    return struct.unpack(">q", struct.pack(">d", x))[0]


def ulp_distance(a, b):
    """Compute ULP distance between two float64 values."""
    if math.isnan(a) or math.isnan(b):
        return math.inf

    if a == b:
        return 0

    ia = float_to_int(a)
    ib = float_to_int(b)

    # Handle sign ordering
    if ia < 0:
        ia = 0x8000000000000000 - ia
    if ib < 0:
        ib = 0x8000000000000000 - ib

    return abs(ia - ib)


# ============================================================
# Characterization tests
# ============================================================


@pytest.mark.parametrize(
    "a,b",
    [
        (0.0, 0.0),
        (1.0, 2.0),
        (-1.0, -2.0),
        (100.0, 100.0),
        (-100.0, -100.0),
        (100.0, -100.0),
        (-1000.0, -1000.0),
    ],
)
def test_matches_reference(a, b):
    expected = reference_lse2(a, b)
    result = robosample.robo_bindings.calculate_log_sum_exp2(a, b)
    assert math.isclose(result, expected, rel_tol=1e-12, abs_tol=1e-12)


def test_symmetry():
    for a, b in [(3.5, -2.1), (100, -100), (-1e3, -1e2)]:
        assert robosample.robo_bindings.calculate_log_sum_exp2(
            a, b
        ) == robosample.robo_bindings.calculate_log_sum_exp2(b, a)


def test_monotonicity():
    a = 1.0
    assert robosample.robo_bindings.calculate_log_sum_exp2(
        a, 3.0
    ) > robosample.robo_bindings.calculate_log_sum_exp2(a, 2.0)


def test_bounds():
    for a, b in [(1, 2), (-10, -20), (100, 0)]:
        res = robosample.robo_bindings.calculate_log_sum_exp2(a, b)
        m = max(a, b)
        assert m <= res <= m + math.log(2)


# ============================================================
# Edge cases
# ============================================================


def test_negative_infinity():
    assert (
        robosample.robo_bindings.calculate_log_sum_exp2(-math.inf, -math.inf)
        == -math.inf
    )


def test_one_infinite():
    assert robosample.robo_bindings.calculate_log_sum_exp2(0.0, -math.inf) == 0.0
    assert robosample.robo_bindings.calculate_log_sum_exp2(-math.inf, 5.0) == 5.0


def test_positive_infinity():
    assert robosample.robo_bindings.calculate_log_sum_exp2(math.inf, 1.0) == math.inf


# ============================================================
# Naive comparison (sanity)
# ============================================================


def naive_lse(a, b):
    return math.log(math.exp(a) + math.exp(b))


@pytest.mark.parametrize(
    "a,b",
    [
        (1000, 1000),
        (-1000, -1000),
        (1000, -1000),
    ],
)
def test_stability_vs_naive(a, b):
    stable = robosample.robo_bindings.calculate_log_sum_exp2(a, b)

    try:
        naive = naive_lse(a, b)
    except (OverflowError, ValueError):
        return

    if math.isfinite(naive):
        assert math.isclose(stable, naive, rel_tol=1e-12)


# ============================================================
# Stress test
# ============================================================


def test_random_stress():
    rng = np.random.default_rng(42)
    for _ in range(10000):
        a = rng.uniform(-1000, 1000)
        b = rng.uniform(-1000, 1000)

        res = robosample.robo_bindings.calculate_log_sum_exp2(a, b)
        ref = reference_lse2(a, b, prec=80)

        assert math.isfinite(res)
        assert math.isclose(res, ref, rel_tol=1e-10, abs_tol=1e-10)


# ============================================================
# Hypothesis strategies
# ============================================================

finite_floats = st.floats(
    min_value=-1e308,
    max_value=1e308,
    allow_nan=False,
    allow_infinity=False,
)

finite_or_inf = st.one_of(
    finite_floats,
    st.just(float("inf")),
    st.just(float("-inf")),
)

# ============================================================
# Property-based tests
# ============================================================


@settings(max_examples=5000)
@given(a=finite_or_inf, b=finite_or_inf)
def test_matches_high_precision(a, b):
    result = robosample.robo_bindings.calculate_log_sum_exp2(a, b)
    ref = reference_lse2(a, b)

    if math.isinf(ref):
        assert result == ref
    else:
        assert math.isfinite(result)
        assert math.isclose(result, ref, rel_tol=1e-10, abs_tol=1e-10)


@settings(max_examples=3000)
@given(a=finite_or_inf, b=finite_or_inf)
def test_symmetry_property(a, b):
    assert robosample.robo_bindings.calculate_log_sum_exp2(
        a, b
    ) == robosample.robo_bindings.calculate_log_sum_exp2(b, a)


@settings(max_examples=3000)
@given(a=finite_or_inf, b=finite_or_inf, c=st.floats(-1e6, 1e6))
def test_translation_invariance(a, b, c):
    if any(math.isinf(x) for x in (a, b, c)):
        return

    lhs = robosample.robo_bindings.calculate_log_sum_exp2(a + c, b + c)
    rhs = robosample.robo_bindings.calculate_log_sum_exp2(a, b) + c

    assert math.isclose(lhs, rhs, rel_tol=1e-10, abs_tol=1e-10)


@settings(max_examples=3000)
@given(a=finite_or_inf, b=finite_or_inf)
def test_bounds_property(a, b):
    res = robosample.robo_bindings.calculate_log_sum_exp2(a, b)
    m = max(a, b)

    if math.isinf(m):
        assert res == m
        return

    assert m <= res <= m + math.log(2)


@settings(max_examples=2000)
@given(
    a=st.floats(0, 1e6),
    delta=st.floats(50, 1e6),
)
def test_dominance(a, delta):
    b = a - delta
    res = robosample.robo_bindings.calculate_log_sum_exp2(a, b)

    assert math.isclose(res, a, rel_tol=1e-12, abs_tol=1e-12)
