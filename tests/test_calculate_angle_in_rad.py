import numpy as np
import robosample
from hypothesis import given, settings
from hypothesis import strategies as st

# -------------------------
# Helpers
# -------------------------


def vecs(min_val=-1e12, max_val=1e12):
    return st.lists(
        st.floats(min_val, max_val, allow_nan=False, allow_infinity=False),
        min_size=3,
        max_size=3,
    ).map(np.array)


# -------------------------
# Deterministic tests
# -------------------------


def test_zero_angle_collinear_same_direction():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    assert np.isclose(angle, 0.0)


def test_pi_angle_collinear_opposite_direction():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([-1.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    assert np.isclose(angle, np.pi)


def test_right_angle():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([0.0, 1.0, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    assert np.isclose(angle, np.pi / 2)


def test_translation_invariance():
    shift = np.array([10.0, -5.0, 3.0])

    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 2.0, 3.0])
    p2 = np.array([-1.0, 4.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    angle2 = robosample.robo_bindings.calculate_angle_in_rad(
        p0 + shift, p1 + shift, p2 + shift
    )

    assert np.isclose(angle1, angle2)


def test_scaling_invariance():
    scale = 7.5

    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 2.0, 3.0])
    p2 = np.array([-1.0, 4.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    angle2 = robosample.robo_bindings.calculate_angle_in_rad(p0, scale * p1, scale * p2)

    assert np.isclose(angle1, angle2)


def test_symmetry():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 2.0, 3.0])
    p2 = np.array([-1.0, 4.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    angle2 = robosample.robo_bindings.calculate_angle_in_rad(p0, p2, p1)

    assert np.isclose(angle1, angle2)


def test_zero_vector_returns_zero():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = p0.copy()
    p2 = np.array([1.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    assert angle == 0.0


def test_both_vectors_zero():
    p0 = np.array([1.0, 1.0, 1.0])
    p1 = p0.copy()
    p2 = p0.copy()

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)
    assert angle == 0.0


# -------------------------
# Property-based tests
# -------------------------


@settings(max_examples=1000)
@given(p0=vecs(), p1=vecs(), p2=vecs())
def test_random_angles_are_valid(p0, p1, p2):
    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)

    assert not np.isnan(angle)
    assert 0.0 <= angle <= np.pi


@settings(max_examples=2000)
@given(
    scale=st.floats(1e-6, 1e12, allow_nan=False, allow_infinity=False),
    eps=st.floats(1e-16, 1e-6, allow_nan=False, allow_infinity=False),
)
def test_near_parallel_vectors(scale, eps):
    p0 = np.zeros(3)
    base = np.array([1.0, 0.0, 0.0]) * scale

    p1 = base
    p2 = base + np.array([0.0, eps * scale, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)

    assert not np.isnan(angle)
    assert 0.0 <= angle <= np.pi


@settings(max_examples=2000)
@given(
    scale=st.floats(1e-6, 1e12, allow_nan=False, allow_infinity=False),
    eps=st.floats(1e-16, 1e-6, allow_nan=False, allow_infinity=False),
)
def test_near_antiparallel_vectors(scale, eps):
    p0 = np.zeros(3)
    base = np.array([1.0, 0.0, 0.0]) * scale

    p1 = base
    p2 = -base + np.array([0.0, eps * scale, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)

    assert not np.isnan(angle)
    assert 0.0 <= angle <= np.pi


@settings(max_examples=1000)
@given(
    scale=st.floats(0.0, 1e-12, allow_nan=False, allow_infinity=False),
    direction=vecs(-1.0, 1.0),
)
def test_near_zero_vectors(scale, direction):
    p0 = np.zeros(3)
    p1 = scale * direction
    p2 = np.array([1.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)

    assert not np.isnan(angle)
    assert 0.0 <= angle <= np.pi


# -------------------------
# Large magnitude (single check sufficient)
# -------------------------


def test_large_magnitude_vectors():
    scale = 1e12

    p0 = np.zeros(3)
    p1 = np.array([1.0, 2.0, 3.0]) * scale
    p2 = np.array([-4.0, 5.0, -6.0]) * scale

    angle = robosample.robo_bindings.calculate_angle_in_rad(p0, p1, p2)

    assert not np.isnan(angle)
    assert 0.0 <= angle <= np.pi
