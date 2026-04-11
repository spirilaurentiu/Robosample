"""
Tests for align_flip_and_translate_frame_along_x_axis
======================================================

What this function does
-----------------------
Given:
  - G_X_F1 : a (3, 4) numpy array representing frame F1 in frame G.
              The first 3 columns are the rotation matrix R (columns = F1's axes in G),
              the last column is the translation p (F1's origin in G).
  - G_v1   : a (3,) numpy array, a point ("station") expressed in G.

The function builds a new frame F3 such that (all expressed in F1's coordinates):

  1. F3's origin  =  v1 expressed in F1
  2. F3's X-axis  =  unit vector pointing FROM v1 BACK TO F1's origin ("flipped")
  3. F3's Y-axis  =  aligned with F1's X-axis (achieved via a dihedral rotation around F3's X-axis)

The return value is F1_X_F3 - i.e. F3's pose expressed in F1's frame - as a (3, 4) array.

Naming convention (SimTK style):
  X_AB  means "transform that expresses frame B as seen from frame A"
  R_AB  is the 3x3 rotation block (columns = B's axes in A)
  p_AB  is the 3-vector translation (B's origin measured from A's origin, expressed in A)

Coordinate frames:
  G   = global/reference frame (the "world")
  F1  = input frame (atom / body frame), specified by G_X_F1
  F3  = output frame, attached at v1, pointing back toward F1

"""

import math

import numpy as np
import pytest
import robosample

# ===========================================================================
# Helpers
# ===========================================================================


def make_transform(R: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Pack a 3x3 rotation R and translation p into a (3,4) array."""
    T = np.zeros((3, 4), dtype=float)
    T[:3, :3] = R
    T[:3, 3] = p
    return T


def identity_transform(p=None) -> np.ndarray:
    """Identity rotation with optional translation."""
    return make_transform(np.eye(3), np.zeros(3) if p is None else np.asarray(p, float))


def rotation_z(angle_rad: float) -> np.ndarray:
    """3x3 rotation matrix around the Z-axis."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def rotation_y(angle_rad: float) -> np.ndarray:
    """3x3 rotation matrix around the Y-axis."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[c, 0.0, s], [0.0, 1.0, 0.0], [-s, 0.0, c]])


def rotation_x(angle_rad: float) -> np.ndarray:
    """3x3 rotation matrix around the X-axis."""
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    return np.array([[1.0, 0.0, 0.0], [0.0, c, -s], [0.0, s, c]])


def extract_R(T: np.ndarray) -> np.ndarray:
    """Extract the 3x3 rotation block from a (3,4) transform."""
    return T[:3, :3]


def extract_p(T: np.ndarray) -> np.ndarray:
    """Extract the translation vector from a (3,4) transform."""
    return T[:3, 3]


def is_valid_rotation(R: np.ndarray, tol: float = 1e-6) -> bool:
    """Return True iff R is an orthogonal matrix with determinant +1."""
    ortho_err = np.max(np.abs(R @ R.T - np.eye(3)))
    det_err = abs(np.linalg.det(R) - 1.0)
    return ortho_err < tol and det_err < tol


def v1_in_F1(G_X_F1: np.ndarray, G_v1: np.ndarray) -> np.ndarray:
    """
    Express G_v1 in F1's coordinate frame.

    Implementation mirrors the C++ code's first two lines:
      G_F1v1  = G_v1 - G_X_F1.p()
      F1_F1v1 = ~G_X_F1.R() * G_F1v1   (~ = transpose = R^{-1} for rotation)
    """
    R_BF = extract_R(G_X_F1)
    p_BF = extract_p(G_X_F1)
    return R_BF.T @ (G_v1 - p_BF)


def expected_x_axis_in_F1(G_X_F1: np.ndarray, G_v1: np.ndarray) -> np.ndarray:
    """
    The expected X-axis of F3 expressed in F1: a unit vector pointing FROM v1
    TOWARD F1's origin (i.e. opposite the F1→v1 direction).
    """
    F1_v1 = v1_in_F1(G_X_F1, G_v1)
    return -F1_v1 / np.linalg.norm(F1_v1)


# ===========================================================================
# Section 1: Output shape and type sanity
# ===========================================================================


class TestOutputShape:
    """The function must return a (3, 4) float array in all cases."""

    def test_returns_numpy_array(self):
        T = identity_transform()
        v = np.array([0.0, 2.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert isinstance(result, np.ndarray), "Result must be a numpy array"

    def test_output_shape_is_3x4(self):
        T = identity_transform()
        v = np.array([0.0, 2.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert result.shape == (3, 4), f"Expected (3,4), got {result.shape}"

    def test_output_is_float(self):
        T = identity_transform()
        v = np.array([0.0, 2.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.issubdtype(result.dtype, np.floating), (
            f"Expected float dtype, got {result.dtype}"
        )


# ===========================================================================
# Section 2: The rotation block must always be a valid rotation matrix
# ===========================================================================


class TestRotationValidity:
    """
    The 3x3 rotation block of the output must always be orthogonal with det=+1.
    If this fails it means F3 is not a proper right-handed frame.
    """

    @pytest.mark.parametrize(
        "v1",
        [
            np.array([0.0, 2.0, 0.0]),  # v1 along +Y from F1
            np.array([0.0, 0.0, 3.0]),  # v1 along +Z from F1
            np.array([0.0, -1.5, 0.0]),  # v1 along -Y from F1
            np.array([0.0, 1.0, 1.0]),  # v1 at 45° in YZ plane
            np.array([0.0, 1.0, -1.0]),
            np.array([0.5, 1.0, 0.5]),  # general direction (not along ±X)
        ],
    )
    def test_rotation_valid_identity_frame(self, v1):
        """With identity F1, output rotation must be valid for all non-degenerate v1."""
        T = identity_transform()
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v1
        )
        R = extract_R(result)
        assert is_valid_rotation(R), (
            f"Output rotation is not valid for v1={v1}\n  R =\n{R}"
        )

    @pytest.mark.parametrize(
        "angle",
        [
            math.pi / 6,  # 30°
            math.pi / 4,  # 45°
            math.pi / 3,  # 60°
            math.pi / 2,  # 90°
            3 * math.pi / 4,
        ],
    )
    def test_rotation_valid_rotated_F1(self, angle):
        """With a rotated F1 frame, output rotation must still be valid."""
        T = make_transform(rotation_z(angle), np.array([1.0, 2.0, 3.0]))
        v = np.array([4.0, 3.0, 0.5])  # arbitrary point not on F1's X-axis
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        R = extract_R(result)
        assert is_valid_rotation(R), (
            f"Output rotation invalid at F1-rotation angle={math.degrees(angle):.0f}°"
        )

    @pytest.mark.parametrize(
        "tx, ty, tz",
        [
            (0.0, 0.0, 0.0),
            (10.0, 0.0, 0.0),
            (-5.0, 3.0, 7.0),
            (100.0, -100.0, 100.0),
        ],
    )
    def test_rotation_valid_translated_F1(self, tx, ty, tz):
        """F1 at arbitrary positions in G; the rotation must still be valid."""
        T = identity_transform(p=[tx, ty, tz])
        # v1 displaced from F1's origin in the Y direction (safe, non-degenerate)
        v = np.array([tx, ty + 2.0, tz])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        R = extract_R(result)
        assert is_valid_rotation(R), (
            f"Output rotation invalid for F1 at ({tx},{ty},{tz})"
        )


# ===========================================================================
# Section 3: Origin of F3 must equal v1 expressed in F1
# ===========================================================================


class TestOriginPlacement:
    """
    The translation part of the output (F3's origin expressed in F1) must equal
    the vector from F1's origin to v1, expressed in F1's own coordinate axes.

    In other words: F1_X_F3.p() == R_F1G * (G_v1 - G_X_F1.p())
    """

    @pytest.mark.parametrize(
        "G_v1",
        [
            np.array([0.0, 3.0, 0.0]),
            np.array([0.0, 0.0, 5.0]),
            np.array([0.0, -2.0, -2.0]),
            np.array([0.0, 1.5, 0.5]),
        ],
    )
    def test_origin_identity_frame(self, G_v1):
        """With identity F1, F3's origin in F1 must equal G_v1 exactly."""
        T = identity_transform()
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        p_out = extract_p(result)
        expected_p = v1_in_F1(T, G_v1)
        np.testing.assert_allclose(
            p_out,
            expected_p,
            atol=1e-10,
            err_msg=f"Origin wrong for G_v1={G_v1}; got {p_out}, expected {expected_p}",
        )

    def test_origin_with_rotated_and_translated_F1(self):
        """
        F1 rotated 45° around Z and placed at [1, 2, 3].
        F3's origin in F1 must be v1 expressed in F1's rotated axes.
        """
        R = rotation_z(math.pi / 4)
        p = np.array([1.0, 2.0, 3.0])
        T = make_transform(R, p)
        G_v1 = np.array([2.0, 3.5, 3.0])

        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        p_out = extract_p(result)
        expected_p = v1_in_F1(T, G_v1)

        np.testing.assert_allclose(
            p_out,
            expected_p,
            atol=1e-10,
            err_msg="Origin wrong for rotated+translated F1",
        )

    def test_origin_distance_preserved(self):
        """
        The distance from F1's origin to v1 must equal the norm of F3's origin in F1.
        This is a scale-invariant property.
        """
        T = make_transform(rotation_y(1.1), np.array([0.0, 0.5, -1.0]))
        G_v1 = np.array([2.0, -1.0, 3.0])

        true_dist = np.linalg.norm(G_v1 - extract_p(T))  # distance in G
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        out_dist = np.linalg.norm(extract_p(result))  # norm of F3 origin in F1

        np.testing.assert_allclose(
            out_dist,
            true_dist,
            rtol=1e-10,
            err_msg="Distance from F1 origin to v1 not preserved",
        )


# ===========================================================================
# Section 4: F3's X-axis must point from v1 TOWARD F1 (the "flip")
# ===========================================================================


class TestXAxisFlip:
    """
    The X-axis of F3 (column 0 of the output rotation) must point in the
    direction from v1 back toward F1's origin.

    In F1 coordinates this direction is  -normalize(F1_F1v1).
    The dot product between the output X-axis and this expected direction
    must equal +1 (they are parallel and unit).
    """

    @pytest.mark.parametrize(
        "G_v1, description",
        [
            (np.array([0.0, 2.0, 0.0]), "v1 along +Y from F1"),
            (np.array([0.0, 0.0, 3.0]), "v1 along +Z from F1"),
            (np.array([0.0, -1.0, 0.0]), "v1 along -Y from F1"),
            (np.array([0.0, 1.0, 1.0]), "v1 in YZ 45°"),
            (np.array([0.0, 1.0, -1.0]), "v1 in YZ -45°"),
            (np.array([0.5, 1.0, 0.5]), "v1 in general direction"),
        ],
    )
    def test_x_axis_points_toward_F1(self, G_v1, description):
        """
        F3.X must be anti-parallel to the (F1 → v1) vector.
        We check alignment (dot product = 1) rather than exact components,
        which is more robust to dihedral rotations around X.
        """
        T = identity_transform()
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )

        x_out = extract_R(result)[:, 0]  # column 0 = F3's X-axis in F1
        x_expected = expected_x_axis_in_F1(T, G_v1)

        alignment = np.dot(x_out, x_expected)
        assert abs(alignment - 1.0) < 1e-6, (
            f"[{description}] F3's X-axis not pointing toward F1.\n"
            f"  got      : {x_out}\n"
            f"  expected : {x_expected}\n"
            f"  dot      : {alignment:.8f} (should be 1.0)"
        )

    def test_x_axis_unit_length(self):
        """F3's X-axis must have unit length."""
        T = identity_transform()
        v = np.array([0.0, 1.5, 2.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        x_col = extract_R(result)[:, 0]
        np.testing.assert_allclose(np.linalg.norm(x_col), 1.0, atol=1e-10)

    def test_x_axis_with_rotated_frame(self):
        """
        Non-identity F1: F3's X-axis must still anti-align with F1→v1, even
        when F1 has a non-trivial orientation and position in G.
        """
        R = rotation_y(math.pi / 3)  # F1 rotated 60° around Y
        T = make_transform(R, np.array([-1.0, 2.0, 0.5]))
        G_v1 = np.array([3.0, 4.0, 1.0])

        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        x_out = extract_R(result)[:, 0]
        x_exp = expected_x_axis_in_F1(T, G_v1)

        np.testing.assert_allclose(
            np.dot(x_out, x_exp),
            1.0,
            atol=1e-6,
            err_msg="X-axis flip incorrect for rotated F1",
        )


# ===========================================================================
# Section 5: F3's Y-axis alignment with F1's X-axis
# ===========================================================================


class TestYAxisAlignment:
    """
    After rotating around F3's X-axis by the computed dihedral angle, F3's
    Y-axis should be aligned with F1's X-axis as much as geometrically possible
    (i.e. the projection of F3.Y onto the plane perpendicular to F3.X should
    maximally align with F1.X).

    Concretely: the dihedral angle between (F1.X, F3.X) measured around F3.X
    should be zero after the transform is applied.
    """

    def _dihedral_around_bond(self, a, b, c, d):
        """
        Compute the dihedral angle of the quadruple a–b–c–d.
        All vectors are 3-D points; the bond axis is b→c.
        Returns angle in radians in [-pi, pi].
        """
        b1 = b - a
        b2 = c - b
        b3 = d - c
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        if np.linalg.norm(n1) < 1e-12 or np.linalg.norm(n2) < 1e-12:
            return 0.0
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        return math.atan2(y, x)

    @pytest.mark.parametrize(
        "G_v1, description",
        [
            (np.array([0.0, 2.0, 0.0]), "v1 along +Y"),
            (np.array([0.0, 0.0, 3.0]), "v1 along +Z"),
            (np.array([0.0, 1.0, 1.0]), "v1 diagonal YZ"),
            (np.array([0.5, 2.0, -0.5]), "v1 general"),
        ],
    )
    def test_y_axis_aligned_with_F1_x(self, G_v1, description):
        """
        F3's Y-axis (in F1 coords) must lie in the plane spanned by F3's X-axis
        and F1's X-axis (i.e. its component perpendicular to F3.X should point
        along F1.X).
        """
        T = identity_transform()
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )

        R_out = extract_R(result)
        F3_x = R_out[:, 0]  # F3's X in F1
        F3_y = R_out[:, 1]  # F3's Y in F1
        F1_x = np.array([1.0, 0.0, 0.0])  # F1's X in F1

        # Project F3.Y and F1.X onto the plane perp to F3.X
        def perp(v, axis):
            return v - np.dot(v, axis) * axis

        F3_y_perp = perp(F3_y, F3_x)
        F1_x_perp = perp(F1_x, F3_x)

        if np.linalg.norm(F1_x_perp) < 1e-9:
            # F1's X is parallel to F3's X: alignment is trivially undefined; skip
            pytest.skip("F1.X parallel to F3.X - Y alignment undefined")

        # They should be parallel (dot product of unit vectors = ±1; "aligned" means +1)
        F3_y_perp_hat = F3_y_perp / np.linalg.norm(F3_y_perp)
        F1_x_perp_hat = F1_x_perp / np.linalg.norm(F1_x_perp)
        dot = np.dot(F3_y_perp_hat, F1_x_perp_hat)

        assert abs(dot - 1.0) < 1e-6, (
            f"[{description}] F3.Y is not aligned with F1.X.\n"
            f"  dot of perp components: {dot:.8f} (expected 1.0)\n"
            f"  F3.X={F3_x}, F3.Y={F3_y}, F1.X={F1_x}"
        )


# ===========================================================================
# Section 6: Known closed-form characterization cases
# ===========================================================================


class TestKnownCases:
    """
    Fully worked-out examples for which we can compute the expected output by hand.
    These pin the exact numerical behaviour and will catch subtle sign / axis / angle bugs.
    """

    def test_identity_frame_v1_along_y(self):
        """
        F1 = identity at origin.  v1 = [0, d, 0].

        Hand-derived:
          F1_F1v1 = [0, d, 0]
          rotAngle  = pi/2  (angle between [0,1,0] and [1,0,0])
          rotAxis   = cross([0,1,0], [1,0,0]) = [0, 0, -1]
          final_angle = -pi/2 + pi = pi/2
          R_F2      = Rot(pi/2, [0,0,-1]) = [[0,1,0],[-1,0,0],[0,0,1]]

          After dihedral: rotates F2 by 0 around its X → F3 = F2.

        Expected output (F1_X_F3):
          origin = [0, 2, 0]
          X-axis (col 0) = [0, -1, 0]  ← from v1 back to F1
          Y-axis (col 1) = [1,  0, 0]  ← aligned with F1's X
          Z-axis (col 2) = [0,  0, 1]  ← right-handed completion
        """
        d = 2.0
        T = identity_transform()
        v = np.array([0.0, d, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )

        R_out = extract_R(result)
        p_out = extract_p(result)

        np.testing.assert_allclose(
            p_out, [0.0, d, 0.0], atol=1e-10, err_msg="Origin wrong for v1 along +Y"
        )
        np.testing.assert_allclose(
            R_out[:, 0],
            [0.0, -1.0, 0.0],
            atol=1e-6,
            err_msg="X-axis wrong for v1 along +Y",
        )
        np.testing.assert_allclose(
            R_out[:, 1],
            [1.0, 0.0, 0.0],
            atol=1e-6,
            err_msg="Y-axis wrong for v1 along +Y",
        )

    def test_identity_frame_v1_along_z(self):
        """
        F1 = identity at origin.  v1 = [0, 0, d].

        Hand-derived:
          F1_F1v1 = [0, 0, d]
          rotAngle  = pi/2
          rotAxis   = cross([0,0,1], [1,0,0]) = [0, 1, 0]
          final_angle = pi/2
          R_F2 col 0 (X-axis) = [0, 0, -1]  ← from v1 back to F1

        Expected:
          origin = [0, 0, d]
          X-axis (col 0) = [0, 0, -1]
        """
        d = 3.0
        T = identity_transform()
        v = np.array([0.0, 0.0, d])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )

        p_out = extract_p(result)
        x_axis = extract_R(result)[:, 0]

        np.testing.assert_allclose(
            p_out, [0.0, 0.0, d], atol=1e-10, err_msg="Origin wrong for v1 along +Z"
        )
        np.testing.assert_allclose(
            x_axis, [0.0, 0.0, -1.0], atol=1e-6, err_msg="X-axis wrong for v1 along +Z"
        )

    def test_identity_frame_v1_along_neg_y(self):
        """
        F1 = identity at origin.  v1 = [0, -d, 0].
        X-axis of F3 must point in +Y (from v1=[0,-d,0] back to F1=[0,0,0]).
        """
        d = 1.5
        T = identity_transform()
        v = np.array([0.0, -d, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )

        p_out = extract_p(result)
        x_axis = extract_R(result)[:, 0]

        np.testing.assert_allclose(
            p_out, [0.0, -d, 0.0], atol=1e-10, err_msg="Origin wrong for v1 along -Y"
        )
        np.testing.assert_allclose(
            x_axis,
            [0.0, 1.0, 0.0],
            atol=1e-6,
            err_msg="X-axis wrong for v1 along -Y (should point back to F1)",
        )

    def test_translated_F1_same_relative_geometry(self):
        """
        Translating F1 and v1 by the same offset must give the same rotation
        and a proportionally shifted origin.

        If G_X_F1 has origin at p0 and v1 = p0 + delta, then shifting both
        by offset changes only the absolute G-coordinates; the relative geometry
        (and thus the rotation block) must be identical.
        """
        delta = np.array([0.0, 2.0, 0.0])
        offset = np.array([10.0, -5.0, 3.0])

        T0 = identity_transform()
        v0 = delta.copy()
        R0 = extract_R(
            robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(T0, v0)
        )

        T1 = identity_transform(p=offset)
        v1 = offset + delta
        R1 = extract_R(
            robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(T1, v1)
        )

        np.testing.assert_allclose(
            R0,
            R1,
            atol=1e-10,
            err_msg="Rotation changed when F1 and v1 were equally translated",
        )

    def test_scaling_distance_only_changes_origin_magnitude(self):
        """
        Scaling v1 farther from F1 (same direction) should only scale F3's
        origin - the rotation must remain identical.
        """
        T = identity_transform()
        v_near = np.array([0.0, 1.0, 0.0])
        v_far = np.array([0.0, 5.0, 0.0])  # same direction, 5x farther

        R_near = extract_R(
            robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                T, v_near
            )
        )
        R_far = extract_R(
            robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                T, v_far
            )
        )

        np.testing.assert_allclose(
            R_near,
            R_far,
            atol=1e-8,
            err_msg="Rotation changed when distance was scaled (same direction)",
        )

        p_near = np.linalg.norm(
            extract_p(
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v_near
                )
            )
        )
        p_far = np.linalg.norm(
            extract_p(
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v_far
                )
            )
        )
        np.testing.assert_allclose(
            p_far,
            5.0 * p_near,
            rtol=1e-10,
            err_msg="Origin magnitude not scaled proportionally",
        )


# ===========================================================================
# Section 7: Consistency / invertibility properties
# ===========================================================================


class TestConsistencyProperties:
    """
    Higher-level properties that must hold regardless of specific input values.
    """

    def test_x_axis_is_anti_parallel_to_F1_origin_direction(self):
        """
        F3.X must be exactly anti-parallel to the vector (F1_origin - v1) expressed
        in F1, which is the same as being parallel to -(v1 - F1_origin) = -F1_F1v1.
        Test with a general, arbitrary configuration.
        """
        R = rotation_z(0.7) @ rotation_y(-0.4)  # compound rotation
        T = make_transform(R, np.array([2.0, -1.0, 0.5]))
        G_v1 = np.array([3.5, 1.0, -0.5])

        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        x_out = extract_R(result)[:, 0]
        x_exp = expected_x_axis_in_F1(T, G_v1)

        np.testing.assert_allclose(
            x_out,
            x_exp,
            atol=1e-6,
            err_msg="X-axis not anti-parallel to F1→v1 for arbitrary config",
        )

    def test_axes_are_mutually_orthogonal(self):
        """X, Y, Z columns of the output rotation must be mutually orthogonal."""
        T = make_transform(rotation_y(1.2), np.array([0.5, -0.5, 1.5]))
        G_v1 = np.array([1.0, 2.0, -1.0])

        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        R = extract_R(result)

        np.testing.assert_allclose(
            np.dot(R[:, 0], R[:, 1]),
            0.0,
            atol=1e-10,
            err_msg="X and Y axes not orthogonal",
        )
        np.testing.assert_allclose(
            np.dot(R[:, 0], R[:, 2]),
            0.0,
            atol=1e-10,
            err_msg="X and Z axes not orthogonal",
        )
        np.testing.assert_allclose(
            np.dot(R[:, 1], R[:, 2]),
            0.0,
            atol=1e-10,
            err_msg="Y and Z axes not orthogonal",
        )

    def test_right_handedness(self):
        """Z-axis must equal X cross Y (right-handed frame)."""
        T = make_transform(rotation_x(0.6), np.zeros(3))
        G_v1 = np.array([0.0, 1.5, -0.5])

        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )
        R = extract_R(result)
        z_expected = np.cross(R[:, 0], R[:, 1])
        np.testing.assert_allclose(
            R[:, 2],
            z_expected,
            atol=1e-10,
            err_msg="Frame is not right-handed (Z ≠ X x Y)",
        )

    def test_output_independent_of_G_frame_rotation(self):
        """
        Rotating both G_X_F1 and G_v1 by the same global rotation R_extra must
        not change the output, because both inputs are just re-expressed in a
        differently-oriented G - the relative geometry is unchanged.

        The output F1_X_F3 is in F1 coordinates, not G coordinates, so applying
        an extra rotation to G should be transparent.
        """
        R_F1 = rotation_z(math.pi / 4)
        p_F1 = np.array([1.0, 0.0, 0.5])
        T = make_transform(R_F1, p_F1)
        G_v1 = np.array([2.0, 2.0, 0.0])
        result0 = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, G_v1
        )

        # Apply an extra rotation of G (rotate everything in G)
        R_extra = rotation_y(math.pi / 3)
        T_rot = make_transform(R_extra @ R_F1, R_extra @ p_F1)
        G_v1_rot = R_extra @ G_v1
        result1 = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T_rot, G_v1_rot
        )

        np.testing.assert_allclose(
            result0,
            result1,
            atol=1e-9,
            err_msg="Output changed when G was rotated - should be invariant",
        )


# ===========================================================================
# Section 8: Stress tests with random inputs
# ===========================================================================


class TestStress:
    """
    Run the function on many random valid inputs and check that all structural
    properties hold. These tests are probabilistic but cover the parameter space
    broadly.
    """

    N_CASES = 200

    def _random_rotation(self, rng: np.random.Generator) -> np.ndarray:
        """Random rotation matrix via QR decomposition of a random matrix."""
        H = rng.standard_normal((3, 3))
        Q, _ = np.linalg.qr(H)
        if np.linalg.det(Q) < 0:
            Q[:, 0] *= -1
        return Q

    def _random_nondegenerate_v1(
        self, rng: np.random.Generator, R_F1: np.ndarray, p_F1: np.ndarray
    ) -> np.ndarray:
        """
        Generate G_v1 that is NOT along F1's ±X axis and NOT coincident with F1's origin.

        Degenerate directions:
          - v1 == F1_origin                     → zero vector crash
          - (v1 - F1_origin) ∥ F1's X-axis     → cross product ≈ 0, degenerate axis
        We avoid these by construction.
        """
        F1_x = R_F1[:, 0]  # F1's X-axis in G
        for _ in range(100):
            # Random displacement in F1 coordinates with zero X-component
            F1_disp = np.array([0.0, rng.uniform(-3.0, 3.0), rng.uniform(-3.0, 3.0)])
            if np.linalg.norm(F1_disp[1:]) < 0.1:
                continue  # too close to degenerate
            G_disp = R_F1 @ F1_disp
            # Add a small random X component to keep it interesting but non-degenerate
            G_disp += F1_x * rng.uniform(-0.5, 0.5)
            v1_candidate = p_F1 + G_disp
            F1_v1 = R_F1.T @ (v1_candidate - p_F1)
            cos_angle = F1_v1[0] / np.linalg.norm(F1_v1)
            # Avoid F1v1 within 10° of ±X in F1 (degenerate cross product)
            if abs(cos_angle) < math.cos(math.radians(10)):
                return v1_candidate
        raise RuntimeError("Could not generate non-degenerate v1 in 100 attempts")

    def test_rotation_always_valid(self):
        """Rotation block must be valid for 200 random inputs."""
        rng = np.random.default_rng(seed=42)
        failures = []
        for i in range(self.N_CASES):
            R = self._random_rotation(rng)
            p = rng.uniform(-5.0, 5.0, 3)
            T = make_transform(R, p)
            v = self._random_nondegenerate_v1(rng, R, p)
            result = (
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v
                )
            )
            Rout = extract_R(result)
            if not is_valid_rotation(Rout):
                failures.append(i)
        assert not failures, (
            f"Rotation invalid in {len(failures)}/{self.N_CASES} random cases"
        )

    def test_origin_always_correct(self):
        """Origin of output must equal v1 in F1 for 200 random inputs."""
        rng = np.random.default_rng(seed=123)
        for i in range(self.N_CASES):
            R = self._random_rotation(rng)
            p = rng.uniform(-5.0, 5.0, 3)
            T = make_transform(R, p)
            v = self._random_nondegenerate_v1(rng, R, p)
            result = (
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v
                )
            )
            p_out = extract_p(result)
            p_exp = v1_in_F1(T, v)
            if not np.allclose(p_out, p_exp, atol=1e-8):
                raise AssertionError(
                    f"Case {i}: origin mismatch.\n"
                    f"  got      : {p_out}\n"
                    f"  expected : {p_exp}"
                )

    def test_x_axis_always_flipped(self):
        """F3's X-axis must point from v1 toward F1 for 200 random inputs."""
        rng = np.random.default_rng(seed=256)
        for i in range(self.N_CASES):
            R = self._random_rotation(rng)
            p = rng.uniform(-5.0, 5.0, 3)
            T = make_transform(R, p)
            v = self._random_nondegenerate_v1(rng, R, p)
            result = (
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v
                )
            )
            x_out = extract_R(result)[:, 0]
            x_exp = expected_x_axis_in_F1(T, v)
            dot = np.dot(x_out, x_exp)
            if abs(dot - 1.0) > 1e-5:
                raise AssertionError(
                    f"Case {i}: X-axis not pointing toward F1.\n"
                    f"  dot={dot:.8f} (expected 1.0)\n"
                    f"  x_out={x_out}, x_exp={x_exp}"
                )

    def test_right_handed_always(self):
        """Output frame must always be right-handed (Z = X x Y)."""
        rng = np.random.default_rng(seed=789)
        for i in range(self.N_CASES):
            R = self._random_rotation(rng)
            p = rng.uniform(-5.0, 5.0, 3)
            T = make_transform(R, p)
            v = self._random_nondegenerate_v1(rng, R, p)
            result = (
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v
                )
            )
            Rout = extract_R(result)
            det = np.linalg.det(Rout)
            if abs(det - 1.0) > 1e-6:
                raise AssertionError(
                    f"Case {i}: det(R) = {det:.8f} - frame not right-handed"
                )

    def test_no_nan_or_inf(self):
        """Output must never contain NaN or Inf for well-conditioned inputs."""
        rng = np.random.default_rng(seed=999)
        for i in range(self.N_CASES):
            R = self._random_rotation(rng)
            p = rng.uniform(-5.0, 5.0, 3)
            T = make_transform(R, p)
            v = self._random_nondegenerate_v1(rng, R, p)
            result = (
                robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
                    T, v
                )
            )
            assert np.all(np.isfinite(result)), (
                f"Case {i}: output contains NaN or Inf\n{result}"
            )


# ===========================================================================
# Section 9: Edge cases - where the function is expected to struggle or fail
# ===========================================================================


class TestEdgeCases:
    """
    These tests document the known degenerate inputs.

    The function uses SimTK::UnitVec3 to normalise the cross product that
    defines the rotation axis.  When that cross product is zero, normalisation
    will produce NaN or crash.

    Two geometrically degenerate situations are:
      A) v1 == F1's origin                    → zero displacement vector
      B) v1 lies exactly on F1's +X ray       → cross([1,0,0],[1,0,0]) = 0

    A third, less obvious case:
      C) v1 lies exactly on F1's −X ray       → cross([-1,0,0],[1,0,0]) = 0
         BUT here the rotation angle = 0, so SimTK may (or may not) handle it.

    The tests below mark these as xfail.  If your PI adds guards that handle
    the degenerate cases gracefully, flip the marks to normal assertions.
    """

    @pytest.mark.xfail(
        strict=False,
        reason="v1 coincides with F1's origin: zero-length vector, "
        "UnitVec normalisation will NaN or crash.",
    )
    def test_v1_equals_F1_origin_crashes(self):
        """
        When v1 == F1's origin, G_F1v1 = 0.  The function attempts to normalise
        it via UnitVec3, which divides by zero.
        Expected behaviour (desired): raise an exception or return NaN.
        Observed behaviour (actual): likely crash or silent NaN.
        """
        T = identity_transform()
        v = np.array([0.0, 0.0, 0.0])  # same as F1's origin
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        # If we reach here, the function returned something - it should be NaN.
        assert not np.all(np.isfinite(result)), (
            "Expected NaN/Inf when v1 == F1 origin, but got finite output"
        )

    @pytest.mark.xfail(
        strict=False,
        reason="v1 exactly on F1's +X axis: cross([1,0,0],[1,0,0])=0, "
        "rotation axis is undefined (angle = pi around zero vector).",
    )
    def test_v1_on_F1_positive_x_axis_crashes(self):
        """
        When v1 is directly along F1's +X axis, the cross product used to find
        the rotation axis is the zero vector.  The function then attempts to
        normalise it as a UnitVec3 (undefined).

        This is the most dangerous degenerate case because the rotation angle
        is pi (not 0), so SimTK cannot silently ignore the axis.

        Example: identity F1, v1 = [2, 0, 0].
        """
        T = identity_transform()
        v = np.array([2.0, 0.0, 0.0])  # exactly on +X of F1
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert not np.all(np.isfinite(result)), (
            "Expected NaN/Inf when v1 is on F1's +X axis"
        )

    def test_v1_near_F1_positive_x_axis_is_fragile(self):
        """
        v1 just 1° off F1's +X axis: the rotation axis exists but the cross
        product is very small (~0.017), so the normalised axis may be inaccurate.

        This is NOT expected to crash, but numerical quality may be poor.
        The test only checks that the output is finite and the rotation is valid.
        """
        T = identity_transform()
        eps = math.radians(1.0)
        v = np.array([math.cos(eps) * 2.0, math.sin(eps) * 2.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.all(np.isfinite(result)), (
            "Output contained NaN/Inf for v1 just 1° off +X axis"
        )
        # Rotation block may be slightly off; just check it's roughly valid
        R = extract_R(result)
        assert is_valid_rotation(R, tol=1e-4), (
            "Rotation not valid for v1 near +X axis (1°)"
        )

    @pytest.mark.xfail(
        strict=False,
        reason="v1 exactly on F1's -X axis: cross([-1,0,0],[1,0,0])=0 "
        "(zero rotation axis). SimTK may handle angle=0 gracefully, "
        "or may not - behaviour is implementation-defined.",
    )
    def test_v1_on_F1_negative_x_axis(self):
        """
        When v1 is directly along F1's −X axis, the rotation angle = (−π + π) = 0.
        The rotation axis is still zero, but a zero rotation is identity, so SimTK
        might handle it silently.  We flag this as xfail because it is geometrically
        degenerate regardless.
        """
        T = identity_transform()
        v = np.array([-2.0, 0.0, 0.0])  # exactly on -X of F1
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        # If it doesn't crash: origin should be at [-2, 0, 0], X-axis at [+1, 0, 0]
        np.testing.assert_allclose(
            extract_p(result),
            [-2.0, 0.0, 0.0],
            atol=1e-6,
            err_msg="Origin wrong for v1 on -X axis",
        )
        x_axis = extract_R(result)[:, 0]
        np.testing.assert_allclose(
            x_axis,
            [1.0, 0.0, 0.0],
            atol=1e-6,
            err_msg="X-axis should point back along +X for v1 on -X axis",
        )

    def test_v1_near_F1_negative_x_axis_works(self):
        """
        v1 just 1° off F1's −X axis: angle ≈ 0, axis nearly zero.
        SimTK's Rotation(small_angle, axis) is usually robust here.
        """
        T = identity_transform()
        eps = math.radians(1.0)
        v = np.array([-math.cos(eps) * 2.0, math.sin(eps) * 2.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.all(np.isfinite(result)), "NaN/Inf for v1 near -X axis"
        R = extract_R(result)
        assert is_valid_rotation(R, tol=1e-4), "Rotation not valid for v1 near -X axis"

    def test_very_small_distance(self):
        """
        v1 very close to but not at F1's origin.
        The function may have precision issues but should not crash.
        """
        T = identity_transform()
        v = np.array([0.0, 1e-8, 0.0])  # 10 nm - tiny but well-defined direction
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.all(np.isfinite(result)), "NaN/Inf for very small F1→v1 distance"

    def test_very_large_distance(self):
        """
        v1 very far from F1 - should not cause overflow.
        """
        T = identity_transform()
        v = np.array([0.0, 1e6, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.all(np.isfinite(result)), "NaN/Inf for very large F1→v1 distance"

    def test_v1_on_F1_y_axis_is_fine(self):
        """Sanity: v1 along F1's Y-axis is the canonical non-degenerate case."""
        T = identity_transform()
        v = np.array([0.0, 3.0, 0.0])
        result = robosample.robo_bindings.align_flip_and_translate_frame_along_x_axis(
            T, v
        )
        assert np.all(np.isfinite(result)), "Unexpected NaN/Inf for safe input"
        assert is_valid_rotation(extract_R(result)), "Rotation invalid for safe input"
