import robosample.robo_bindings as rb
import pytest
import numpy as np


# ============================================================
# 1. CORE CHIRALITY TESTS (triple product orientation)
# ============================================================

def test_chirality_same():
    """
    Identical source and target frames -> same chirality -> no mismatch.
    """
    s1, s2, s3 = np.eye(3)
    assert not rb.is_chirality_mismatch(s1,s2,s3,s1,s2,s3)


def test_chirality_inverted():
    """
    Reflection (flip one axis) -> determinant sign changes -> mismatch.
    """
    s1, s2, s3 = np.eye(3)
    t3 = -s3
    assert rb.is_chirality_mismatch(s1,s2,s3,s1,s2,t3)


def test_reflection_flip():
    """
    Full reflection matrix should invert chirality.
    """
    s = np.eye(3)
    t = np.diag([1,1,-1]) @ s
    assert rb.is_chirality_mismatch(s[0],s[1],s[2], t[0],t[1],t[2])


def test_random_chirality_consistency():
    """
    Random frames copied exactly -> must never detect mismatch.
    Guards against sign instability.
    """
    for _ in range(1000):
        s = np.random.randn(3,3)
        assert not rb.is_chirality_mismatch(s[0],s[1],s[2],
                                            s[0],s[1],s[2])


# ============================================================
# 2. NUMERICAL EDGE CASES (stability)
# ============================================================

def test_non_unit_vectors():
    """
    Chirality depends only on orientation, not magnitude.
    Scaling vectors must not affect result.
    """
    s1 = np.array([2,0,0])
    s2 = np.array([0,3,0])
    s3 = np.array([0,0,4])
    assert not rb.is_chirality_mismatch(s1,s2,s3,s1,s2,s3)


def test_chirality_near_zero():
    """
    Nearly collinear vectors -> triple product ~ 0.
    Should NOT trigger mismatch due to numerical noise.
    """
    eps = 1e-10
    s1 = np.array([1,0,0])
    s2 = np.array([1,eps,0])
    s3 = np.array([0,0,1])
    assert not rb.is_chirality_mismatch(s1,s2,s3,s1,s2,s3)


# ============================================================
# 3. PLANARITY / GEOMETRY TESTS
# ============================================================

def test_planarity_not_broken():
    """
    Vector lies in plane -> deviation small -> below threshold.
    """
    v = np.array([1,0,0])
    n = np.array([0,0,1])
    dev = rb.signed_plane_deviation(v, n)
    assert not rb.exceeds_planarity_threshold(dev, 0.01)


def test_planarity_broken():
    """
    Vector aligned with normal -> maximal deviation -> must break.
    """
    v = np.array([0,0,1])
    n = np.array([0,0,1])
    dev = rb.signed_plane_deviation(v, n)
    assert rb.exceeds_planarity_threshold(dev, 0.01)


def test_threshold_sensitivity():
    """
    Construct vector with exact angle θ:
        v = (cosθ, 0, sinθ)
    so deviation = sinθ exactly.
    """
    theta = 0.01
    v = np.array([np.cos(theta), 0, np.sin(theta)])
    n = np.array([0,0,1])

    dev = rb.signed_plane_deviation(v, n)

    assert not rb.exceeds_planarity_threshold(dev, 0.1)
    assert rb.exceeds_planarity_threshold(dev, 0.001)


# ============================================================
# 4. REFERENCE INDEX RESOLUTION (topological mapping)
# ============================================================

def test_reference_indices_basic():
    """
    Identity mapping: BC indices already ordered.
    """
    idx = rb.resolve_reference_indices([0,1,2])
    assert idx.zero == 0
    assert idx.one == 1
    assert idx.two == 2


def test_reference_indices_permuted():
    """
    Permuted BC indices -> must correctly locate 0 and 1.
    """
    idx = rb.resolve_reference_indices([2,0,1])
    assert idx.zero == 1
    assert idx.one == 2


def test_reference_indices_duplicate():
    """
    Duplicate BC indices -> should not crash.
    Behavior: assign first valid occurrence.
    """
    idx = rb.resolve_reference_indices([0,0,1])
    assert idx.zero in [0,1]
    assert idx.one == 2


def test_reference_indices_missing():
    """
    Missing BC 0/1 -> fallback to defaults.
    """
    idx = rb.resolve_reference_indices([2,3,4])
    assert idx.zero == 0
    assert idx.one == 1


# ============================================================
# 5. PER-BOND CHIRALITY (local comparison)
# ============================================================

def test_per_bond_mismatch():
    """
    Local inversion of a single bond -> mismatch detected.
    """
    s1, s2 = [1,0,0], [0,1,0]
    si = [0,0,1]
    ti = [0,0,-1]

    assert rb.is_bond_chirality_mismatch(s1,s2,si,s1,s2,ti)


def test_per_bond_degenerate():
    """
    Degenerate geometry -> should not produce false mismatch.
    """
    eps = 1e-12
    s1 = [1,0,0]
    s2 = [1,eps,0]
    si = [0,0,1]

    assert not rb.is_bond_chirality_mismatch(s1,s2,si,s1,s2,si)


# ============================================================
# 6. PLANAR -> CHIRAL MAPPING
# ============================================================

def test_plane_positive():
    """
    Positive deviation -> RightHanded.
    """
    assert rb.chirality_from_plane_deviation(0.1) == rb.BondCenterChirality.RightHanded


def test_plane_negative():
    """
    Negative deviation -> LeftHanded.
    """
    assert rb.chirality_from_plane_deviation(-0.1) == rb.BondCenterChirality.LeftHanded


def test_plane_zero():
    """
    Zero deviation -> currently mapped to RightHanded (design choice).
    """
    assert rb.chirality_from_plane_deviation(0.0) == rb.BondCenterChirality.RightHanded


# ============================================================
# 7. CHIRALITY FLIPPING
# ============================================================

def test_flip_right():
    """Right -> Left"""
    assert rb.flipped_chirality(rb.BondCenterChirality.RightHanded) == rb.BondCenterChirality.LeftHanded


def test_flip_left():
    """Left -> Right"""
    assert rb.flipped_chirality(rb.BondCenterChirality.LeftHanded) == rb.BondCenterChirality.RightHanded


def test_flip_planar():
    """Planar remains unchanged"""
    assert rb.flipped_chirality(rb.BondCenterChirality.Planar) == rb.BondCenterChirality.Planar


# ============================================================
# 8. INTEGRATION (minimal pipeline sanity)
# ============================================================

def test_full_chirality_pipeline():
    """
    Simple 3-vector inversion -> pipeline detects mismatch.
    """
    s = np.eye(3)
    t = np.diag([1,1,-1]) @ s

    assert rb.is_chirality_mismatch(s[0],s[1],s[2],
                                    t[0],t[1],t[2])
    
# ============================================================
# 9. GEOMETRIC INVARIANCE & NUMERICAL BOUNDARIES
# ============================================================

def test_mixed_scaling_invariance():
    """
    Chirality must be invariant under arbitrary scaling of vectors.

    Even if each axis is scaled differently (including sign flips),
    chirality depends only on orientation (sign of triple product),
    not magnitude.

    Here:
    - s = canonical basis
    - t = scaled + two sign flips -> net positive orientation

    Expected: NO mismatch.
    """
    s1 = np.array([1,0,0])
    s2 = np.array([0,2,0])
    s3 = np.array([0,0,3])

    t1 = 10*s1
    t2 = -5*s2
    t3 = -2*s3  # two sign flips -> preserve chirality

    assert not rb.is_chirality_mismatch(s1,s2,s3,t1,t2,t3)


def test_rotation_invariance():
    """
    Chirality must be invariant under proper rotations (determinant = +1).

    Applying a rotation matrix should not change orientation sign.

    This guards against implementations that accidentally depend on
    coordinate frame instead of relative geometry.
    """
    R = np.array([
        [0,-1,0],
        [1, 0,0],
        [0, 0,1]
    ])  # 90° rotation about z-axis

    s = np.eye(3)
    t = R @ s

    assert not rb.is_chirality_mismatch(
        s[0],s[1],s[2],
        t[0],t[1],t[2]
    )


def test_exact_coplanar():
    """
    Degenerate case: vectors lie in the same plane -> triple product = 0.

    Chirality is undefined in this case, but implementation should:
    - NOT report mismatch
    - remain numerically stable

    Guards against division-by-zero or noisy sign flips.
    """
    s1 = [1,0,0]
    s2 = [0,1,0]
    s3 = [1,1,0]  # coplanar

    assert not rb.is_chirality_mismatch(s1,s2,s3,s1,s2,s3)


def test_threshold_boundary():
    """
    Boundary condition: deviation exactly at threshold.

    Tests consistency of comparison:
        abs(dev) >= sin(threshold)

    Ensures no off-by-one / floating-point inconsistency
    at the decision boundary.
    """
    n = np.array([0,0,1])
    v = np.array([0,0,np.sin(0.01)])

    dev = rb.signed_plane_deviation(v, n)

    assert rb.exceeds_planarity_threshold(dev, 0.01)