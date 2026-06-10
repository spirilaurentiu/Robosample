#!/usr/bin/env python3
"""
kinematic_soa.py  --  SoA O(n) Rigid Body Kinematics   [v3]
=============================================================

Algorithm  : Featherstone RBDA (position + velocity FK)
Complexity : O(n_bodies) traversal | O(n_atoms) atom scatter (vectorized)
Layout     : Structure-of-Arrays -- flat numpy arrays, zero per-body objects
Joints     : Weld | Pin | Prismatic | Free  (int enum, no vtable)
H-matrix   : motion subspace, computed procedurally from (joint_type, axis)
Traversal  : BFS topological sort, computed once at finalize()
Visualizer : pyvista (preferred) or matplotlib fallback

------------------------------------------------------------------------------
JAIN vs FEATHERSTONE -- May 2026
------------------------------------------------------------------------------
Both are O(n) and equivalent in complexity.

Jain's Spatial Operator Algebra (JPL 1991, book 2011) frames everything as
tip-to-base operator products (phi operators, inertia operators, factored
Newton-Euler).  Elegant for deriving new algorithms and proving correctness.

Featherstone's RBDA notation maps more directly to SIMD/GPU code.  Pinocchio 3
(LAAS-CNRS, de-facto reference in 2026), Drake, and MuJoCo all use Featherstone
recursion internally.

For 100k rigid bodies / 1M atoms:
  Bodies are likely many INDEPENDENT small trees (CG molecules), not one
  100k-DOF serial chain.  The serial O(n) recursion per tree is trivially fast;
  the bottleneck is throughput across trees -- batched FK.
  cuRobo (NVIDIA 2023+) and Isaac Lab batch FK/dynamics across GPU threads.
  On CPU: numba.prange across trees, serial O(n) within each tree.

Verdict: Featherstone for implementation, Jain for reasoning about the
Jacobian-transpose force projection Q = J^T F.

------------------------------------------------------------------------------
SIMBODY COMPARISON (MatterSubsystem)
------------------------------------------------------------------------------
Simbody models:
    qdot = N u                       kinematic differential equations
    M udot + G^T mult = f            equations of motion
    G udot = b                       constraint equations

What Simbody has that we DO NOT (and why):

  N matrix (qdot = N u):
    For PIN and PRISMATIC joints, N = I so qdot = u directly.
    For quaternion FREE joints, N is 4x3.  We currently assume N = I,
    which is wrong for proper quaternion integration of free bodies.
    Flag: FREE joint q should eventually be a 7-vector (quaternion + xyz).

  Constraint subsystem G = [P; V; A]:
    Holonomic (position), non-holonomic (velocity), acceleration constraints.
    We have no constraints -- not needed for our MD use case.

  System Jacobian J (partial velocity matrix, explicit):
    We compute J^T F implicitly via the Phi tip-to-base pass instead.

  Composite body inertia R_i:
    Tip-to-base sum assuming all joints below are locked.
    Needed for recursive Newton-Euler (inverse dynamics).
    NOT needed for ABA (forward dynamics) -- we skip it.

  Articulated body inertia I_A[i]:
    Tip-to-base pass of ABA; accounts for joint compliance.
    Needed for O(n) forward dynamics.  Stubbed here.

  Quaternion normalization constraints n(q):
    Needed when FREE joints use quaternions.  Not yet implemented.

What Simbody has that we DO have (equivalent):
  - Tree topology (parent array, BFS topological order)
  - Position FK: X_GB via base-to-tip transform accumulation
  - Velocity FK: V_GB via Phi^T base-to-tip recursion (this file)
  - Phi matrix concept (computed inline from body origins)
  - Per-body mass
  - Joint types and H-matrix

What we have that Simbody does not:
  - SoA layout (Simbody is AoS + PIMPL)
  - Atom stations and O(n_atoms) scatter
  - Direct OpenMM force pipeline hooks

------------------------------------------------------------------------------
REALIZE STAGES -- what they are and whether we need them
------------------------------------------------------------------------------
Simbody stages and our equivalents:

  realizeInstance     : allocation / one-time setup
                        -> our __init__ + finalize()
                        STATUS: done.

  realizePosition     : q -> X_GB, Phi, COM_G, Mk_G
                        -> our forward_kinematics()
                        STATUS: done (X_GB + atom scatter).
                        MISSING: stored Phi array, full Mk_G spatial inertia.

  realizeVelocity     : (q, u) -> V_GB, coriolis a[i], gyroscopic b[i]
                        -> our realize_velocity()
                        STATUS: V_GB implemented.
                        MISSING: coriolis and gyroscopic (deferred to ABA).
                        NEEDED FOR: kinetic energy, integrator, ABA bias.

  realizeDynamics     : articulated body inertia pass (tip -> base, ABA Pass 2)
                        -> our _realize_articulated_body_inertia() [STUB]
                        STATUS: not implemented.
                        NEEDED FOR: O(n) forward dynamics.

  realizeAcceleration : ABA Pass 3 (base -> tip), computes udot and A_GB
                        -> our _realize_acceleration() [STUB]
                        STATUS: not implemented.
                        NEEDED FOR: forward dynamics.

  Composite body inertia (separate from ABA):
                        Tip-to-base sum of locked-joint inertias.
                        Needed for Newton-Euler inverse dynamics.
                        NOT needed for our ABA + OpenMM path -- skipped.

  Priority order for future implementation:
    1. realize_velocity: done here.
    2. ABA pass 2 + 3: when adding an integrator.
    3. N matrix fix for FREE joints: when quaternion integration matters.
    4. Full Mk_G: when kinetic energy or gyroscopic forces are needed.
"""

# ---- imports -----------------------------------------------------------------
import time
from enum import IntEnum
from typing import List, Optional

import numpy as np

import robosample

# ---- optional numba ----------------------------------------------------------
try:
    import numba as nb
    from numba import njit, prange

    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

# ---- visualization: pyvista preferred, matplotlib fallback -------------------
try:
    import pyvista as pv

    HAS_PYVISTA = True
except ImportError:
    HAS_PYVISTA = False
    try:
        import matplotlib

        matplotlib.use("TkAgg")
    except Exception:
        pass
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# ==============================================================================
# TRANSFORM NOTATION
# ==============================================================================
#
# Convention: X_AB is the pose of frame B expressed in frame A.
#   It converts coordinates FROM frame B TO frame A:
#       p_A = X_AB * p_B
#
# A 4x4 homogeneous transform is stored as:
#   X_AB = [ R_AB  t_AB ]
#           [  0     1  ]
#
#   R_AB (3x3 rotation):
#       Columns are the B-frame unit axes expressed in A coordinates.
#       Example: R_GB[:,0] = body X-axis direction expressed in world space.
#       If you have a vector v_body in body coordinates,
#       then v_world = R_GB @ v_body.
#
#   t_AB (3x1 translation):
#       Position of B's origin expressed in A coordinates.
#       Example: t_GB = where the body origin sits in world space.
#
# Full transform:
#   p_world = R_GB @ p_body + t_GB
#   (this is exactly what the atom scatter einsum computes)
#
# How we store it:
#   X_world[i]          = X_GB for body i, shape (4,4)
#   X_world[i, :3, :3]  = R_GB  (rotation)
#   X_world[i, :3,  3]  = t_GB  (translation, body origin in world)
#
# ------------------------------------------------------------------------------
# FRAME NAMES (following Simbody convention, simplified)
# ------------------------------------------------------------------------------
#
#   G = Ground  (world frame, inertial, fixed to the simulation box)
#   P = Parent body frame
#   B = Body (child) frame
#   F = Parent joint frame -- fixed to P, offset to where the joint lives on P
#   M = Child joint frame  -- fixed to B, offset to where the joint lives on B
#
# Full Simbody chain:
#   G -> P -> F -> [joint DOF] -> M -> B
#   X_PB = X_PF * X_FM(q) * X_MB    (three frames, one moving piece X_FM)
#
# What we do (COLLAPSED):
#   G -> P -> [joint DOF] -> B
#   T_local  = joint_transform(q) @ X_ref
#   X_world[i] = X_world[parent] @ T_local
#   where X_ref ~ X_PF * X_MB  (the two fixed offsets merged into one)
#
# Consequence of collapsing: we cannot separately recover the joint frame
# positions (F and M).  Fine for atom scatter and force projection.
# If joint-frame Jacobians are needed later, split X_ref into X_PF + X_MB.
#
# ------------------------------------------------------------------------------
# PHI MATRIX (shift operator)
# ------------------------------------------------------------------------------
#
#   Phi(r) = [ I   0 ]   (6x6 spatial matrix)
#             [ r*  I ]
#   where r* = skew(r) = cross-product matrix of the 3-vector r.
#   r = t_world[child] - t_world[parent]  (in Ground frame)
#
#   Phi^T shifts spatial velocities OUTWARD (parent -> child):
#       V_child = Phi^T * V_parent + H * u
#
#   Phi   shifts spatial forces   INWARD  (child -> parent):
#       F_parent += Phi * F_child
#
#   We compute Phi inline (no explicit 6x6 storage needed).
# ==============================================================================


# ==============================================================================
# JOINT MODEL
#
# Spatial vector layout: [omega_x, omega_y, omega_z, v_x, v_y, v_z]
#   Indices [0:3] = angular (rotation) part
#   Indices [3:6] = linear  (velocity) part
#   (Featherstone section 2.2 convention)
#
# H-matrix (Featherstone calls it S; Jain calls it H-hat):  shape (6, nDOF)
#   Encodes which spatial DOFs the joint leaves free.
#   Spatial velocity recursion:
#       V_GB[child] = Phi^T(r) * V_GB[parent] + H * u
# ==============================================================================


class JointType(IntEnum):
    WELD = 0  # 0 DOF -- rigid weld, no motion allowed
    PIN = 1  # 1 DOF -- revolute: rotation about a fixed axis
    PRISMATIC = 2  # 1 DOF -- prismatic: translation along a fixed axis
    FREE = 3  # 6 DOF -- unconstrained floating body
    # CAUTION: q is a scalar (rotation angle) in this prototype.
    # Production use requires q as 7-vector (unit quaternion + xyz).


JOINT_DOF = np.array([0, 1, 1, 6], dtype=np.int32)

# CPK element colors: atomic number -> RGB float [0..1]
CPK = {
    1: np.array([0.90, 0.90, 0.90]),  # H  white
    6: np.array([0.22, 0.22, 0.22]),  # C  dark grey
    7: np.array([0.20, 0.40, 0.85]),  # N  blue
    8: np.array([0.85, 0.15, 0.15]),  # O  red
    15: np.array([1.00, 0.50, 0.00]),  # P  orange
    16: np.array([0.90, 0.82, 0.10]),  # S  yellow
    17: np.array([0.15, 0.75, 0.15]),  # Cl green
}
_DEFAULT_COLOR = np.array([0.60, 0.60, 0.70])


# ---- spatial math ------------------------------------------------------------


def _rodrigues(axis: np.ndarray, theta: float) -> np.ndarray:
    """
    3x3 rotation matrix via Rodrigues formula.
    axis need not be unit-length (normalized internally).
    """
    n = axis / (np.linalg.norm(axis) + 1e-15)
    K = np.array(
        [[0, -n[2], n[1]], [n[2], 0, -n[0]], [-n[1], n[0], 0]], dtype=np.float64
    )
    return np.eye(3) + np.sin(theta) * K + (1.0 - np.cos(theta)) * (K @ K)


def _T(R: np.ndarray, p: np.ndarray) -> np.ndarray:
    """Assemble a 4x4 homogeneous transform from 3x3 rotation R and 3-vector p."""
    T = np.eye(4, dtype=np.float64)
    T[:3, :3] = R
    T[:3, 3] = p
    return T


def _Ttrans(p: np.ndarray) -> np.ndarray:
    """Pure-translation 4x4 homogeneous transform (rotation = identity)."""
    T = np.eye(4, dtype=np.float64)
    T[:3, 3] = p
    return T


def joint_transform(jtype: int, axis: np.ndarray, q: float) -> np.ndarray:
    """
    Compute T_joint(q): the 4x4 SE(3) transform contributed by the joint DOF.

    This corresponds to X_FM(q) in Simbody notation -- the moving piece that
    the generalized coordinate q drives.

    The full local transform is:
        T_local = joint_transform(q) @ X_ref
    where X_ref is the constant parent->body offset at q=0 (our collapsed X_PF*X_MB).

    Then:
        X_world[child] = X_world[parent] @ T_local

    PIN:       rotates about `axis` by angle q (Rodrigues formula)
    PRISMATIC: translates along `axis` by distance q
    FREE:      rotates about `axis` by angle q (simplified -- see class note)
    WELD:      identity (no DOF, no motion)
    """
    if jtype == JointType.WELD:
        return np.eye(4, dtype=np.float64)
    elif jtype == JointType.PIN:
        return _T(_rodrigues(axis, q), np.zeros(3))
    elif jtype == JointType.PRISMATIC:
        ax = axis / (np.linalg.norm(axis) + 1e-15)
        return _Ttrans(ax * q)
    elif jtype == JointType.FREE:
        return _T(_rodrigues(axis, q), np.zeros(3))
    return np.eye(4, dtype=np.float64)


def H_matrix(jtype: int, axis: np.ndarray) -> np.ndarray:
    """
    Motion subspace H: shape (6, nDOF).

    Encodes which spatial DOFs the joint leaves free.
    Featherstone calls this S; Jain calls it H-hat.

    Spatial velocity recursion:
        V_GB[child] = Phi^T(r) * V_GB[parent] + H * u

    Layout of H:
        Rows [0:3] correspond to angular (omega) spatial DOFs.
        Rows [3:6] correspond to linear (v) spatial DOFs.

    PIN (revolute about axis):
        H[:,0] = [ax, 0, 0, 0]  -- only the angular rows are non-zero.
        Meaning: the joint contributes angular velocity but no linear velocity
        at the joint frame origin.

    PRISMATIC (translate along axis):
        H[:,0] = [0, 0, 0, ax]  -- only the linear rows are non-zero.

    WELD:
        H is (6, 0) -- zero columns, no free DOF.

    FREE:
        H = I_6  -- all 6 DOFs free.
        (Simplified because q is scalar here; production needs H columns
        corresponding to the quaternion parametrization.)
    """
    ax = axis / (np.linalg.norm(axis) + 1e-15)
    if jtype == JointType.WELD:
        return np.zeros((6, 0), dtype=np.float64)
    elif jtype == JointType.PIN:
        H = np.zeros((6, 1), dtype=np.float64)
        H[:3, 0] = ax  # angular rows = axis
        return H
    elif jtype == JointType.PRISMATIC:
        H = np.zeros((6, 1), dtype=np.float64)
        H[3:, 0] = ax  # linear rows = axis
        return H
    elif jtype == JointType.FREE:
        # PROTOTYPE SIMPLIFICATION: q and u are scalars (rotation about axis).
        # H is therefore (6,1) selecting that angular DOF -- same structure as PIN.
        # Production: q = 7-vector (quaternion + xyz), H = I_6, u = 6D spatial vel.
        H = np.zeros((6, 1), dtype=np.float64)
        H[:3, 0] = ax
        return H
    return np.zeros((6, 0), dtype=np.float64)


# ==============================================================================
# RigidBodySystem -- pure SoA, no per-body heap objects
#
# Body arrays, indexed [0 .. n_bodies):
#   parent         int32  (N,)      index of parent body; -1 = root (child of Ground)
#   jtype          int32  (N,)      JointType enum value
#   jaxis          f64    (N,3)     joint axis in parent frame
#   X_ref          f64    (N,4,4)   constant parent->body offset at q=0 (X_PF*X_MB)
#   q              f64    (N,)      generalized position  (1 scalar per body in prototype)
#   u              f64    (N,)      generalized velocity  (1 scalar per body in prototype)
#   mass           f64    (N,)
#
# Computed by forward_kinematics (realize_position):
#   X_world        f64    (N,4,4)   X_GB: body frame pose in Ground
#                                   X_world[i,:3,:3] = R_GB (rotation)
#                                   X_world[i,:3, 3] = t_GB (origin in world)
#
# Computed by realize_velocity:
#   V_GB           f64    (N,6)     spatial velocity [omega_x,y,z, v_x,y,z] in Ground
#
# Atom arrays, indexed [0 .. n_atoms):
#   ab             int32  (A,)      atom -> body index
#   station        f64    (A,3)     atom position in body frame ("atom station")
#   amass          f64    (A,)
#   element        int32  (A,)      atomic number
#
# Computed by forward_kinematics (atom scatter):
#   pos_world      f64    (A,3)     atom position in Ground frame
#
# Force pipeline (filled externally, consumed by project_forces_to_generalized):
#   F_atom_world   f64    (A,3)     per-atom force in Ground frame (from OpenMM)
#   F_body_spatial f64    (N,6)     [torque; force] per body, Ground frame
#   Q              f64    (N,)      generalized force per joint DOF (J^T F)
# ==============================================================================


class RigidBodySystem:
    """
    Structure-of-Arrays rigid body / atom system.

    Usage:
        rbs = RigidBodySystem()
        b0 = rbs.add_body("base",   parent_id=-1, joint=JointType.WELD, ...)
        b1 = rbs.add_body("link_1", parent_id=b0, joint=JointType.PIN,  ...)
        rbs.add_atoms(b1, stations=np.random.randn(20, 3))
        rbs.finalize()

        # Each MD step:
        rbs.forward_kinematics()       # q -> X_world, pos_world
        rbs.realize_velocity()         # u -> V_GB
        mock_openmm_forces(rbs)        # placeholder for OpenMM call
        scatter_atom_forces(rbs)       # F_atom_world -> F_body_spatial
        project_forces_to_generalized(rbs)  # F_body_spatial -> Q
        # (future) rbs._realize_articulated_body_inertia()  # ABA pass 2
        # (future) rbs._realize_acceleration()              # ABA pass 3 -> udot
    """

    def __init__(self, max_bodies: int = 200_000, max_atoms: int = 2_000_000):
        N = max_bodies
        A = max_atoms

        # -- Body SoA ----------------------------------------------------------
        self.n_bodies = 0
        self.parent = np.full(N, -1, dtype=np.int32)
        self.jtype = np.zeros(N, dtype=np.int32)
        self.jaxis = np.zeros((N, 3), dtype=np.float64)
        self.X_ref = np.broadcast_to(np.eye(4), (N, 4, 4)).copy()
        self.q = np.zeros(N, dtype=np.float64)
        self.u = np.zeros(N, dtype=np.float64)
        self.mass = np.ones(N, dtype=np.float64)
        self.bname = np.empty(N, dtype=object)
        self.bname[:] = ""

        # root_welded: if True, root bodies (parent == -1) are fixed to Ground.
        # If False, the root's joint DOF is live (floating base).
        self.root_welded: bool = True

        # -- Computed by forward_kinematics (realize_position) -----------------
        # X_world[i] = X_GB for body i.
        # X_world[i, :3, :3] = R_GB: rotation matrix, columns = body axes in world.
        # X_world[i, :3,  3] = t_GB: body origin position in world.
        self.X_world = np.broadcast_to(np.eye(4), (N, 4, 4)).copy()

        # -- Computed by realize_velocity --------------------------------------
        # V_GB[i] = 6D spatial velocity of body i in Ground frame.
        # V_GB[i, 0:3] = omega: angular velocity of the body, in Ground
        # V_GB[i, 3:6] = v:     linear velocity of the body ORIGIN, in Ground
        self.V_GB = np.zeros((N, 6), dtype=np.float64)

        # -- Atom SoA ----------------------------------------------------------
        self.n_atoms = 0
        self.ab = np.zeros(A, dtype=np.int32)  # atom -> body
        self.station = np.zeros((A, 3), dtype=np.float64)  # pos in body frame
        self.amass = np.ones(A, dtype=np.float64)
        self.element = np.full(A, 6, dtype=np.int32)  # atomic number

        # -- Computed by forward_kinematics (atom scatter) ---------------------
        # pos_world[j] = R_GB[body[j]] @ station[j] + t_GB[body[j]]
        self.pos_world = np.zeros((A, 3), dtype=np.float64)

        # -- Force pipeline ----------------------------------------------------
        # F_atom_world: per-atom force in Ground frame, written by OpenMM (or mock).
        self.F_atom_world = np.zeros((A, 3), dtype=np.float64)
        # F_body_spatial: accumulated per-body spatial force.
        #   [:, 0:3] = torque about body origin, in Ground
        #   [:, 3:6] = total force on body, in Ground
        self.F_body_spatial = np.zeros((N, 6), dtype=np.float64)
        # Q: generalized force per joint DOF.  Result of J^T F.
        # Feeds into ABA Pass 2 as the applied tau[i].
        self.Q = np.zeros(N, dtype=np.float64)

        # -- Internal ----------------------------------------------------------
        self._topo = np.zeros(N, dtype=np.int32)
        self._ready = False

    # ---- builder API ---------------------------------------------------------

    def add_body(
        self,
        name: str,
        parent_id: int,  # -1 for roots (children of Ground)
        joint: JointType,
        axis: np.ndarray,  # joint axis in parent frame
        X_ref: Optional[np.ndarray] = None,  # 4x4: parent->body at q=0 (X_PF*X_MB)
        mass: float = 1.0,
        q0: float = 0.0,
        u0: float = 0.0,
    ) -> int:
        """Add one body. Returns its body_id (0-indexed integer)."""
        i = self.n_bodies
        self.n_bodies += 1
        self.parent[i] = parent_id
        self.jtype[i] = int(joint)
        self.jaxis[i] = np.asarray(axis, dtype=np.float64)
        if X_ref is not None:
            self.X_ref[i] = np.asarray(X_ref, dtype=np.float64)
        self.mass[i] = mass
        self.q[i] = q0
        self.u[i] = u0
        self.bname[i] = name
        self._ready = False
        return i

    def add_atoms(
        self,
        body_id: int,
        stations: np.ndarray,  # (K, 3) positions in body frame
        masses: Optional[np.ndarray] = None,
        elements: Optional[np.ndarray] = None,
    ) -> slice:
        """Bulk-add atoms to a body. Returns the slice they occupy."""
        K = len(stations)
        j0 = self.n_atoms
        self.n_atoms += K
        self.ab[j0 : j0 + K] = body_id
        self.station[j0 : j0 + K] = stations
        if masses is not None:
            self.amass[j0 : j0 + K] = masses
        if elements is not None:
            self.element[j0 : j0 + K] = elements
        return slice(j0, j0 + K)

    # ---- setup ---------------------------------------------------------------

    def finalize(self):
        """
        BFS topological sort of the body forest.
        Guarantees parents always appear before their children in _topo.
        Call once after all add_body() / add_atoms() calls.
        """
        n = self.n_bodies
        ch = [[] for _ in range(n)]
        q = []
        for i in range(n):
            p = self.parent[i]
            if p < 0:
                q.append(i)
            else:
                ch[p].append(i)
        order, head = [], 0
        while head < len(q):
            node = q[head]
            head += 1
            order.append(node)
            q.extend(ch[node])
        if len(order) != n:
            raise ValueError(
                "Body graph is disconnected -- every body must trace to a root"
            )
        self._topo[:n] = order
        self._ready = True

    # ---- realize_position ----------------------------------------------------

    def forward_kinematics(self):
        """
        Pass 1a (base -> tip): position kinematics.
        Corresponds to Simbody's realizePosition stage.

        For each body i in topological order (parents guaranteed before children):

            T_joint  = joint_transform(q[i], axis[i])   -- X_FM(q): the moving piece
            T_local  = T_joint @ X_ref[i]               -- X_PB: parent->body transform
            X_world[i] = X_world[parent] @ T_local      -- X_GB: world transform

        X_world[i] interpretation:
            X_world[i, :3, :3] = R_GB  -- body orientation in Ground
                                          columns = [body_x, body_y, body_z] in world
            X_world[i, :3,  3] = t_GB  -- body origin position in world

        Atom scatter (vectorized, O(n_atoms)):
            pos_world[j] = R_GB[body[j]] @ station[j] + t_GB[body[j]]
            Implemented as a batched einsum for cache efficiency.
        """
        assert self._ready, "Call finalize() before forward_kinematics()"

        n = self.n_bodies
        na = self.n_atoms

        for idx in range(n):
            i = self._topo[idx]
            p = self.parent[i]

            T_j = joint_transform(self.jtype[i], self.jaxis[i], self.q[i])
            T_local = T_j @ self.X_ref[i]

            if p < 0:
                self.X_world[i] = T_local
            else:
                self.X_world[i] = self.X_world[p] @ T_local

        if na > 0:
            # Bw[j] = 4x4 world transform of atom j's parent body
            Bw = self.X_world[self.ab[:na]]  # (na, 4, 4)
            R = Bw[:, :3, :3]  # (na, 3, 3)
            t = Bw[:, :3, 3]  # (na, 3)
            # pos_world[j] = R[j] @ station[j] + t[j]
            self.pos_world[:na] = np.einsum("bij,bj->bi", R, self.station[:na]) + t

    # ---- realize_velocity ----------------------------------------------------

    def realize_velocity(self):
        """
        Pass 1b (base -> tip): velocity kinematics.
        Corresponds to Simbody's realizeVelocity stage.
        Must be called after forward_kinematics() (needs current X_world).

        Recursion (Featherstone Eq. 5.3):
            V_GB[i] = Phi^T(r_i) * V_GB[parent] + H[i] * u[i]

        Phi^T (outward velocity shift) acts on V_parent = [omega_p, v_p]:
            omega_child = omega_p
            v_child     = v_p - r x omega_p
        where:
            r = t_GB[child] - t_GB[parent]   (vector from parent origin to child origin, in Ground)
            x = cross product

        Physical meaning of the shift:
            If a rigid parent body has angular velocity omega_p and its origin moves
            at v_p, then a point at displacement r from its origin has linear velocity
            v_p + omega_p x r.  Rearranging for the child origin convention gives
            v_child = v_p - r x omega_p.

        H[i] * u[i] -- joint's direct contribution to velocity:
            PIN:       [ax * u, 0, 0, 0]  -- pure angular, no linear at joint origin
            PRISMATIC: [0, 0, 0, ax * u]  -- pure linear, no angular
            WELD:      zero vector

        V_GB[i] storage:
            V_GB[i, 0:3] = omega  (angular velocity in Ground)
            V_GB[i, 3:6] = v      (linear velocity of body ORIGIN in Ground)

        NOT YET COMPUTED (deferred to ABA):
            Coriolis acceleration a[i] = Phi^T * a[parent] + A_mobilizer
            Gyroscopic force b[i] = m * [omega x (J*omega),  omega x (omega x r_com)]
        These are needed for the ABA bias force z[i].
        """
        assert self._ready, "Call finalize() before realize_velocity()"

        n = self.n_bodies

        for idx in range(n):
            i = self._topo[idx]
            p = self.parent[i]

            if p < 0:
                V_parent = np.zeros(6)
            else:
                V_parent = self.V_GB[p]

            omega_p = V_parent[:3]
            v_p = V_parent[3:]

            # r = vector from parent origin to this body's origin, in Ground
            if p < 0:
                r = np.zeros(3)
            else:
                r = self.X_world[i, :3, 3] - self.X_world[p, :3, 3]

            # Phi^T * V_parent
            omega_shifted = omega_p
            v_shifted = v_p - np.cross(r, omega_p)

            # H * u
            H = H_matrix(self.jtype[i], self.jaxis[i])
            V_joint = H @ np.atleast_1d(self.u[i]) if H.shape[1] > 0 else np.zeros(6)

            self.V_GB[i, :3] = omega_shifted + V_joint[:3]
            self.V_GB[i, 3:] = v_shifted + V_joint[3:]

    # ---- ABA stubs (forward dynamics, not yet implemented) ------------------

    def _realize_articulated_body_inertia(self):
        """
        STUB -- Pass 2 (tip -> base): articulated body inertia (ABA).

        Corresponds to Simbody's realizeDynamics stage.
        NOT IMPLEMENTED.  Required for O(n) forward dynamics.

        Algorithm (Featherstone ABA, tip-to-base):
        Start: I_A[i] = Mk_G[i]  (spatial body inertia, 6x6)
        For each body i in REVERSE topological order:
            D[i]     = H[i]^T @ I_A[i] @ H[i]               (dof x dof)
            G[i]     = I_A[i] @ H[i] @ inv(D[i])             (6 x dof)
            P[i]     = I_A[i] - G[i] @ H[i]^T @ I_A[i]      (6x6, projected)
            I_A[p]  += Phi[i] @ P[i] @ Phi[i]^T              (propagate inward)
            z[p]    += Phi[i] @ (z[i] + I_A[i]@c[i] + G[i]@(tau[i] - H[i]^T@z[i]))

        where:
            Mk_G[i]  = spatial inertia of body i in Ground (not yet computed)
            c[i]     = coriolis acceleration (requires realize_velocity)
            z[i]     = velocity-dependent bias force
            tau[i]   = generalized applied force (our Q[i])

        Jain's insight: the product of (I - G[i]*H[i]^T) propagated inward
        is the O(n) factored form of M^{-1}.  This is why ABA is O(n).

        COMPOSITE BODY INERTIA (different thing, NOT needed for ABA):
            R_i = Mk_G[i] + sum_children Phi_j @ R_j @ Phi_j^T
            Assumes all downstream joints are LOCKED (rigid subtree).
            Used for Newton-Euler inverse dynamics.
            We do NOT need inverse dynamics -- skipped permanently.

        Implement when: adding an integrator.
        Requires first: coriolis a[i] and gyroscopic b[i] from realize_velocity,
                        and full Mk_G[i] (spatial inertia) in realize_position.
        """
        raise NotImplementedError(
            "_realize_articulated_body_inertia: ABA Pass 2 not implemented. "
            "See docstring for full algorithm."
        )

    def _realize_acceleration(self):
        """
        STUB -- Pass 3 (base -> tip): acceleration propagation (ABA).

        Corresponds to Simbody's realizeAcceleration stage.
        NOT IMPLEMENTED.  Required for O(n) forward dynamics.

        Algorithm (Featherstone ABA, base-to-tip):
        Start: A_GB[root] = -gravity  (expressed as 6D spatial vector)
        For each body i in topological order:
            udot[i]    = inv(D[i]) @ (tau[i] - H[i]^T @ (I_A[i]@A_GB[parent] + z[i]))
            A_GB[i]    = Phi^T[i] @ A_GB[parent] + H[i]@udot[i] + c[i]

        After this pass, udot[i] is the generalized acceleration.

        Integration (simple Euler, for reference):
            u_new = u + udot * dt
            q_new = q + u_new * dt    (valid when N = I, i.e. not quaternions)

        Implement when: adding an integrator.
        Requires first: _realize_articulated_body_inertia (Pass 2).
        """
        raise NotImplementedError(
            "_realize_acceleration: ABA Pass 3 not implemented. "
            "See docstring for full algorithm."
        )

    # ---- diagnostics ---------------------------------------------------------

    def summary(self):
        print("-" * 60)
        print(f"  RigidBodySystem  |  {self.n_bodies} bodies  |  {self.n_atoms} atoms")
        print(f"  root_welded = {self.root_welded}")
        print("-" * 60)
        print(
            f"  {'idx':>4}  {'name':<14}  {'parent':>6}  {'joint':<10}  {'axis':<22}  q"
        )
        print(f"  {'-' * 4}  {'-' * 14}  {'-' * 6}  {'-' * 10}  {'-' * 22}  ---")
        for i in range(min(self.n_bodies, 20)):
            jname = JointType(self.jtype[i]).name
            ax = self.jaxis[i]
            print(
                f"  {i:>4}  {str(self.bname[i]):<14}  {self.parent[i]:>6}  "
                f"{jname:<10}  [{ax[0]:+.2f} {ax[1]:+.2f} {ax[2]:+.2f}]  {self.q[i]:.3f}"
            )
        if self.n_bodies > 20:
            print(f"  ... {self.n_bodies - 20} more bodies ...")
        print()
        print("  H-matrices for first 6 bodies:")
        for i in range(min(self.n_bodies, 6)):
            H = H_matrix(self.jtype[i], self.jaxis[i])
            ax = self.jaxis[i]
            dof = H.shape[1]
            print(
                f"    body {i} ({self.bname[i]}, {JointType(self.jtype[i]).name}, "
                f"axis=[{ax[0]:.2f},{ax[1]:.2f},{ax[2]:.2f}]) -> "
                f"H shape {H.shape}, DOF={dof}"
            )
            if dof > 0:
                labels = ["wx", "wy", "wz", "vx", "vy", "vz"]
                for row in range(6):
                    vals = "  ".join(f"{H[row, c]:+.1f}" for c in range(dof))
                    print(f"        {labels[row]}: {vals}")
        print("-" * 60)


# ==============================================================================
# OPENMM FORCE PIPELINE
#
# Each MD step, the full pipeline is:
#
#   rbs.forward_kinematics()              -- q -> X_world, pos_world
#   rbs.realize_velocity()               -- u -> V_GB
#   mock_openmm_forces(rbs)              -- [MOCK] pos_world -> F_atom_world
#   scatter_atom_forces(rbs)             -- F_atom_world -> F_body_spatial
#   project_forces_to_generalized(rbs)   -- F_body_spatial -> Q  (J^T F)
#   rbs._realize_articulated_body_inertia()  [STUB] -- ABA pass 2
#   rbs._realize_acceleration()              [STUB] -- ABA pass 3 -> udot
#   integrate: u += udot*dt, q += u*dt
# ==============================================================================


def mock_openmm_forces(rbs: RigidBodySystem) -> None:
    """
    MOCK PLACEHOLDER for the OpenMM force evaluation.

    In production this function will:
      1. Send rbs.pos_world[:rbs.n_atoms] (shape: n_atoms x 3, in Ground frame)
         to OpenMM via context.setPositions().
      2. Call context.getState(getForces=True).
      3. Read back per-atom forces in Ground frame (same coordinate system).
      4. Write them into rbs.F_atom_world[:rbs.n_atoms].

    This mock does nothing: F_atom_world keeps its current value
    (all zeros after __init__).  The function signature and array
    contract are final and will not change.

    Contract:
      Input:  rbs.pos_world[:rbs.n_atoms]    shape (n_atoms, 3)  Ground frame
      Output: rbs.F_atom_world[:rbs.n_atoms] shape (n_atoms, 3)  Ground frame

    Units: match OpenMM convention (kJ/mol/nm or kcal/mol/A depending on system).
    OpenMM forces are already in Ground (world) frame -- no rotation needed.
    """
    pass  # replace with: rbs.F_atom_world[:rbs.n_atoms] = openmm_context.getForces()


def scatter_atom_forces(rbs: RigidBodySystem) -> None:
    """
    Sum per-atom forces (Ground frame) into per-body spatial forces.

    Must be called after mock_openmm_forces() has filled F_atom_world.

    For each atom j belonging to body i:
        torque contribution:  F_body[i, 0:3] += moment_arm[j] x F_atom[j]
        force  contribution:  F_body[i, 3:6] += F_atom[j]

    moment_arm[j] = pos_world[j] - t_GB[body[j]]
        = vector from body origin to atom position, in Ground frame.

    The cross product moment_arm x force gives the torque that the atom's
    force exerts about the body's origin, expressed in Ground frame.

    Result in rbs.F_body_spatial:
        [:, 0:3] = total torque about body origin, Ground frame
        [:, 3:6] = total force on body, Ground frame

    Cost: O(n_atoms), fully vectorized.
    np.add.at is used for scatter (handles multiple atoms per body correctly).
    """
    na = rbs.n_atoms
    rbs.F_body_spatial[: rbs.n_bodies] = 0.0
    if na == 0:
        return

    body_origins = rbs.X_world[rbs.ab[:na], :3, 3]  # (na, 3)  t_GB per atom
    moment_arms = rbs.pos_world[:na] - body_origins  # (na, 3)  r from origin to atom
    torques = np.cross(moment_arms, rbs.F_atom_world[:na])  # (na, 3)

    np.add.at(rbs.F_body_spatial, (rbs.ab[:na], slice(0, 3)), torques)
    np.add.at(rbs.F_body_spatial, (rbs.ab[:na], slice(3, 6)), rbs.F_atom_world[:na])


def project_forces_to_generalized(rbs: RigidBodySystem) -> None:
    """
    Pass 2 tip->base: project spatial body forces to generalized forces.

    Must be called after scatter_atom_forces() has filled F_body_spatial.

    Computes Q = J^T F via the Phi tip-to-base recursion.
    This avoids forming the full system Jacobian J explicitly (O(n^2)).

    For each body i in REVERSE topological order (leaves first):
        Q[i]      = H[i]^T @ F_body_spatial[i]    -- project to free DOF
        F[parent] += Phi(r) @ F[i]                -- shift force inward

    Phi(r) applied to spatial force F = [tau, f]:
        tau_shifted = tau + r x f    (torque shifts with moment arm)
        f_shifted   = f              (force magnitude is preserved)
    where r = t_GB[child] - t_GB[parent]  (in Ground)

    Physical meaning: the child body's force and torque, when transferred
    to the parent's reference point, gains an additional moment r x f.

    Result in rbs.Q[:n_bodies]:
        Q[i] = generalized force for body i's joint DOF.
        Feeds into ABA Pass 2 as the applied tau[i].
        NOT udot -- that comes from ABA Pass 3.

    This implements Q = J^T F in Jain's notation:
        Q = H^T * (accumulated tip-to-base spatial force)
    O(n) via Phi recursion instead of O(n^2) explicit Jacobian.
    """
    n = rbs.n_bodies
    F = rbs.F_body_spatial[:n].copy()  # working copy; do not mutate original

    for idx in range(n - 1, -1, -1):
        i = rbs._topo[idx]
        p = rbs.parent[i]

        H = H_matrix(rbs.jtype[i], rbs.jaxis[i])
        rbs.Q[i] = float((H.T @ F[i])[0]) if H.shape[1] > 0 else 0.0

        if p >= 0:
            r = rbs.X_world[i, :3, 3] - rbs.X_world[p, :3, 3]
            tau_i = F[i, :3]
            force_i = F[i, 3:]
            F[p, :3] += tau_i + np.cross(r, force_i)
            F[p, 3:] += force_i


# ==============================================================================
# DEMO SCENE BUILDERS
# ==============================================================================


def build_robot_arm(
    rbs: RigidBodySystem,
    n_links: int = 6,
    link_length: float = 1.4,
    n_atoms_per_link: int = 15,
) -> List[int]:
    """
    Serial-chain robot arm welded to Ground.

    Frame layout per link:
        body 0:  base, WELD to Ground, X_ref = I (at world origin)
        body k:  PIN about Z, X_ref = T([link_length, 0, 0])
                 -- joint at parent's origin; body frame origin link_length away.

    At q=0 for all joints: arm lies along the +X world axis.
    """
    rng = np.random.default_rng(7)
    ids = []

    base = rbs.add_body(
        "arm_base",
        parent_id=-1,
        joint=JointType.WELD,
        axis=np.array([0.0, 0.0, 1.0]),
        X_ref=np.eye(4),
        mass=2.0,
    )
    ids.append(base)
    rbs.add_atoms(base, rng.uniform(-0.3, 0.3, (5, 3)), elements=rng.integers(6, 8, 5))

    prev = base
    for k in range(1, n_links + 1):
        X_ref = _Ttrans(np.array([link_length, 0.0, 0.0]))
        body = rbs.add_body(
            f"link_{k}",
            parent_id=prev,
            joint=JointType.PIN,
            axis=np.array([0.0, 0.0, 1.0]),
            X_ref=X_ref,
            mass=1.0,
            q0=np.pi / (n_links + 1) * k * 0.5,
        )
        ids.append(body)

        t_vals = rng.uniform(0.05, 0.85, n_atoms_per_link)
        r_vals = rng.uniform(0.05, 0.20, n_atoms_per_link)
        phi_vals = rng.uniform(0, 2 * np.pi, n_atoms_per_link)
        stations = np.column_stack(
            [t_vals * link_length, r_vals * np.cos(phi_vals), r_vals * np.sin(phi_vals)]
        )
        ep = [6, 7, 8, 6, 6, 15]
        rbs.add_atoms(
            body,
            stations,
            elements=np.array([ep[i % len(ep)] for i in range(n_atoms_per_link)]),
        )
        prev = body

    X_ee = _Ttrans(np.array([link_length, 0.0, 0.0]))
    ee = rbs.add_body(
        "end_effector",
        parent_id=prev,
        joint=JointType.PIN,
        axis=np.array([0.0, 1.0, 0.0]),
        X_ref=X_ee,
        mass=0.5,
        q0=0.0,
    )
    ids.append(ee)
    rbs.add_atoms(ee, rng.standard_normal((8, 3)) * 0.15, elements=np.full(8, 8))
    return ids


def build_molecule_clusters(
    rbs: RigidBodySystem,
    n_molecules: int = 40,
    atoms_per_mol: int = 18,
    spread: float = 9.0,
) -> List[int]:
    """Free-floating rigid bodies representing coarse-grained molecule clusters."""
    rng = np.random.default_rng(13)
    ids = []
    ep = [6, 6, 7, 8, 6, 16, 6, 7, 8, 6, 6, 15, 6, 6, 7, 8, 6, 6, 6, 7]

    for k in range(n_molecules):
        pos = rng.standard_normal(3)
        pos = pos / (np.linalg.norm(pos) + 1e-9) * rng.uniform(2, spread)

        body = rbs.add_body(
            f"mol_{k:03d}",
            parent_id=-1,
            joint=JointType.FREE,
            axis=rng.standard_normal(3),
            X_ref=_T(np.eye(3), pos),
            mass=float(atoms_per_mol) * 12.0,
            q0=rng.uniform(0, 2 * np.pi),
        )
        ids.append(body)

        shell = rng.uniform(-1, 1, (atoms_per_mol - 3, 3))
        shell /= np.linalg.norm(shell, axis=1, keepdims=True) + 1e-9
        shell *= rng.uniform(0.2, 0.55, (atoms_per_mol - 3, 1))
        core = rng.standard_normal((3, 3)) * 0.10
        rbs.add_atoms(
            body,
            np.vstack([shell, core]),
            elements=np.array(
                [ep[i % len(ep)] for i in range(atoms_per_mol)], dtype=np.int32
            ),
        )
    return ids


# ==============================================================================
# ANIMATION
# ==============================================================================


def animate_step(
    rbs: RigidBodySystem, robot_ids: List[int], mol_ids: List[int], t: float
) -> None:
    for k, bid in enumerate(robot_ids):
        if rbs.jtype[bid] == JointType.PIN:
            base_angle = np.pi / len(robot_ids) * k
            rbs.q[bid] = base_angle + 0.45 * np.sin(t * 0.8 + k * 0.6)
    for k, bid in enumerate(mol_ids):
        rbs.q[bid] = t * 0.35 + bid * 0.17


# ==============================================================================
# PYVISTA VISUALIZER
# ==============================================================================


def run_pyvista(rbs: RigidBodySystem, robot_ids: List[int], mol_ids: List[int]) -> None:
    rbs.forward_kinematics()
    rbs.realize_velocity()
    na = rbs.n_atoms
    nb_b = rbs.n_bodies

    colors = np.array(
        [CPK.get(int(e), _DEFAULT_COLOR) for e in rbs.element[:na]], dtype=np.float32
    )

    pl = pv.Plotter(window_size=(1400, 800), title="SoA Rigid Body Kinematics v3")
    pl.set_background("#0d1117")
    pl.add_axes(color="gray")

    cloud = pv.PolyData(rbs.pos_world[:na].copy())
    cloud["colors"] = (colors * 255).astype(np.uint8)
    pl.add_mesh(
        cloud, render_points_as_spheres=True, point_size=7, rgb=True, scalars="colors"
    )

    fcloud = pv.PolyData(rbs.X_world[:nb_b, :3, 3].copy())
    pl.add_mesh(fcloud, render_points_as_spheres=True, point_size=10, color="#00ff88")

    skeleton = pv.lines_from_points(
        np.array([rbs.X_world[b, :3, 3] for b in robot_ids])
    )
    pl.add_mesh(skeleton, color="#ff6b35", line_width=3)

    pl.add_text(
        f"{nb_b} bodies | {na} atoms | Featherstone O(n) FK+vel",
        position="upper_left",
        font_size=11,
        color="#aaaaaa",
    )

    print("\n[pyvista] 3D viewer open -- close window to exit\n")
    pl.show(auto_close=False, interactive_update=True)
    t0 = time.time()
    try:
        while True:
            t = time.time() - t0
            animate_step(rbs, robot_ids, mol_ids, t)
            rbs.forward_kinematics()
            rbs.realize_velocity()
            cloud.points = rbs.pos_world[:na].copy()
            fcloud.points = rbs.X_world[:nb_b, :3, 3].copy()
            skeleton.points = np.array([rbs.X_world[b, :3, 3] for b in robot_ids])
            pl.update()
            time.sleep(1.0 / 60)
    except KeyboardInterrupt:
        pass
    finally:
        pl.close()


# ==============================================================================
# MATPLOTLIB FALLBACK
# ==============================================================================


def run_matplotlib(
    rbs: RigidBodySystem, robot_ids: List[int], mol_ids: List[int]
) -> None:
    rbs.forward_kinematics()
    rbs.realize_velocity()
    na = rbs.n_atoms

    fig = plt.figure(figsize=(13, 8), facecolor="#0d1117")
    ax = fig.add_subplot(111, projection="3d", facecolor="#0d1117")
    ax.set_title("SoA Rigid Body Kinematics v3", color="white", fontsize=13, pad=12)
    ax.tick_params(colors="#555555")
    for pane in [ax.xaxis.pane, ax.yaxis.pane, ax.zaxis.pane]:
        pane.fill = False
    ax.grid(True, color="#222222")

    colors = np.array([CPK.get(int(e), _DEFAULT_COLOR) for e in rbs.element[:na]])
    sc = ax.scatter(
        rbs.pos_world[:na, 0],
        rbs.pos_world[:na, 1],
        rbs.pos_world[:na, 2],
        c=colors,
        s=14,
        alpha=0.85,
        depthshade=True,
    )
    (arm_line,) = ax.plot(
        [rbs.X_world[b, 0, 3] for b in robot_ids],
        [rbs.X_world[b, 1, 3] for b in robot_ids],
        [rbs.X_world[b, 2, 3] for b in robot_ids],
        "o-",
        color="#ff6b35",
        lw=2.5,
        ms=6,
    )
    span = 12
    ax.set_xlim(-span, span)
    ax.set_ylim(-span, span)
    ax.set_zlim(-span, span)
    ax.set_xlabel("X", color="#666")
    ax.set_ylabel("Y", color="#666")
    ax.set_zlabel("Z", color="#666")

    t0 = time.time()

    def _update(_frame):
        t = time.time() - t0
        animate_step(rbs, robot_ids, mol_ids, t)
        rbs.forward_kinematics()
        rbs.realize_velocity()
        sc._offsets3d = (
            rbs.pos_world[:na, 0],
            rbs.pos_world[:na, 1],
            rbs.pos_world[:na, 2],
        )
        arm_line.set_data_3d(
            [rbs.X_world[b, 0, 3] for b in robot_ids],
            [rbs.X_world[b, 1, 3] for b in robot_ids],
            [rbs.X_world[b, 2, 3] for b in robot_ids],
        )
        return (sc, arm_line)

    ani = FuncAnimation(fig, _update, interval=16, blit=False, cache_frame_data=False)
    plt.tight_layout()
    print("\n[matplotlib] 3D window open -- close to exit\n")
    plt.show()


# ==============================================================================
# BENCHMARK
# ==============================================================================


def benchmark(rbs: RigidBodySystem, n_iters: int = 400) -> None:
    for _ in range(10):  # warm-up
        rbs.forward_kinematics()
        rbs.realize_velocity()

    t0 = time.perf_counter()
    for _ in range(n_iters):
        rbs.forward_kinematics()
        rbs.realize_velocity()
    dt = (time.perf_counter() - t0) / n_iters

    print()
    print("-" * 60)
    print(f"  BENCHMARK  ({n_iters} iters, FK + velocity)")
    print(f"  {rbs.n_bodies:>8} bodies  ->  {dt * 1e3:7.3f} ms / step")
    print(
        f"  {rbs.n_atoms:>8} atoms   ->  {rbs.n_atoms / dt / 1e6:7.2f} M atom-transforms / sec"
    )
    print("  Body FK: serial O(n).  Atom scatter: vectorized O(n_atoms).")
    print("  For 100k independent molecules: parallelize across trees")
    print("  with numba.prange or a CUDA batched-FK kernel (cuRobo style).")
    print("-" * 60)

    # smoke-test force pipeline (all zeros -- just verifies shapes and code paths)
    mock_openmm_forces(rbs)
    scatter_atom_forces(rbs)
    project_forces_to_generalized(rbs)
    print(f"  Force pipeline smoke: Q[:5] = {rbs.Q[:5]}  (zero -- mock forces)")
    print("-" * 60)
    print()


# ==============================================================================
# MAIN
# ==============================================================================


def main() -> None:

    print()
    print("=" * 60)
    print("  SoA O(n) Rigid Body Kinematics  [v3]")
    print("  Featherstone FK+vel  |  SoA  |  ASCII-only")
    print("=" * 60)

    context = robosample.Context(
        name="ala-dipeptide",
        seed=6000,
        prmtop="examples/ala-dipeptide.prmtop",
        inpcrd="examples/ala-dipeptide.rst7",
        write_freq=100,
        testing=False,
    )

    rbs = RigidBodySystem(max_bodies=10_000, max_atoms=100_000)
    rbs.root_welded = True

    robot_ids = build_robot_arm(rbs, n_links=args.links, n_atoms_per_link=12)
    mol_ids = build_molecule_clusters(
        rbs, n_molecules=args.mols, atoms_per_mol=args.apm, spread=args.spread
    )
    rbs.finalize()
    rbs.summary()
    benchmark(rbs, n_iters=400)

    if not args.no_viz:
        if HAS_PYVISTA:
            run_pyvista(rbs, robot_ids, mol_ids)
        else:
            run_matplotlib(rbs, robot_ids, mol_ids)
    else:
        print("Visualization skipped (--no-viz).")


if __name__ == "__main__":
    main()
