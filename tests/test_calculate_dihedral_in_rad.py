import numpy as np
import robosample
from hypothesis import given, settings
from hypothesis import strategies as st

# -----------------------------
# Helpers
# -----------------------------


def random_rotation_matrix():
    # Simple random rotation via QR decomposition
    A = np.random.randn(3, 3)
    Q, _ = np.linalg.qr(A)
    if np.linalg.det(Q) < 0:
        Q[:, 0] *= -1
    return Q


def vecs():
    return st.lists(
        st.floats(-1e6, 1e6, allow_nan=False, allow_infinity=False),
        min_size=3,
        max_size=3,
    ).map(np.array)


def angular_distance(a, b):
    return min(abs(a - b), 2 * np.pi - abs(a - b))


# -----------------------------
# Deterministic tests
# -----------------------------


def test_planar_dihedral_zero():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    assert angular_distance(angle, 0.0) < 1e-8


def test_dihedral_pi_non_degenerate():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])

    # Flip in opposite direction but keep coplanarity
    p3 = np.array([3.0, 0.0, 0.0])  # fully collinear → still degenerate

    # Better:
    p2 = np.array([1.0, 1.0, 0.0])
    p3 = np.array([2.0, -1.0, 0.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)

    assert angular_distance(angle, np.pi) < 1e-8


def test_dihedral_pi():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])

    # Reflect across plane defined by (p0, p1, p2)
    p3 = np.array([3.0, 0.0, 0.0])

    # Then introduce perpendicular displacement to ensure well-defined planes
    p2 = np.array([2.0, 1.0, 0.0])
    p3 = np.array([3.0, -1.0, 0.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)

    assert angular_distance(angle, np.pi) < 1e-8


def test_chirality_sign_non_degenerate():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 1.0, 0.0])  # break collinearity

    p3 = np.array([3.0, 1.0, 1.0])
    p3_reflected = np.array([3.0, 1.0, -1.0])

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(
        p0, p1, p2, p3_reflected
    )

    assert not np.isclose(angle1, 0.0)
    assert np.sign(angle1) == -np.sign(angle2)


def test_translation_invariance():
    shift = np.array([5.0, -3.0, 2.0])

    p0 = np.array([0.1, 0.2, 0.3])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 1.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(
        p0 + shift, p1 + shift, p2 + shift, p3 + shift
    )

    assert angular_distance(angle1, angle2) < 1e-8


def test_scaling_invariance():
    scale = 10.0

    p0 = np.array([0.1, 0.2, 0.3])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 1.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(
        p0 * scale, p1 * scale, p2 * scale, p3 * scale
    )

    assert angular_distance(angle1, angle2) < 1e-8


def test_collinear_bonds():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 0.0, 0.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)

    assert not np.isnan(angle)
    assert -np.pi <= angle <= np.pi


def test_rotation_invariance():
    R = random_rotation_matrix()

    p0 = np.array([0.1, 0.2, 0.3])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 1.0, 0.5])

    p_rot = [R @ p for p in (p0, p1, p2, p3)]

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(*p_rot)

    assert angular_distance(angle1, angle2) < 1e-8


def test_periodicity():
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([2.0, 1.0, 1.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)

    assert angular_distance(angle, angle + 2 * np.pi) < 1e-8


# -----------------------------
# Property-based tests
# -----------------------------


@settings(max_examples=1000)
@given(p0=vecs(), p1=vecs(), p2=vecs(), p3=vecs())
def test_dihedral_valid_range(p0, p1, p2, p3):
    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    assert not np.isnan(angle)
    assert -np.pi <= angle <= np.pi


@settings(max_examples=2000)
@given(
    eps=st.floats(1e-16, 1e-6, allow_nan=False, allow_infinity=False),
)
def test_near_collinearity(eps):
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, eps, 0.0])
    p3 = np.array([3.0, eps, 0.0])

    angle = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    assert not np.isnan(angle)
    assert -np.pi <= angle <= np.pi


@settings(max_examples=1000)
@given(vec=vecs())
def test_reflection_flips_sign(vec):
    p0 = np.array([0.0, 0.0, 0.0])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = vec

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)

    p3_reflected = np.array([vec[0], vec[1], -vec[2]])
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(
        p0, p1, p2, p3_reflected
    )

    if abs(angle1) > 1e-8:
        assert angular_distance(angle1, -angle2) < 1e-8


@settings(max_examples=1000)
@given(scale=st.floats(1e-12, 1e12, allow_nan=False, allow_infinity=False))
def test_scale_extremes(scale):
    p0 = np.array([0.1, 0.2, 0.3])
    p1 = np.array([1.0, 0.0, 0.0])
    p2 = np.array([2.0, 0.0, 0.0])
    p3 = np.array([3.0, 1.0, 0.5])

    angle1 = robosample.robo_bindings.calculate_dihedral_in_rad(p0, p1, p2, p3)
    angle2 = robosample.robo_bindings.calculate_dihedral_in_rad(
        p0 * scale, p1 * scale, p2 * scale, p3 * scale
    )

    assert angular_distance(angle1, angle2) < 1e-8
