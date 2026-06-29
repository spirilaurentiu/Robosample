# The mobilizer frame chain (Simbody convention)

Each mobilized body `B` is positioned by a chain of frames between ground (`G`) and `B`. Following the Simbody convention, the frames are:

- **Frame_G** -- Ground (inertial) frame.
- **Frame_P** -- parent body frame (the inboard body's frame).
- **Frame_F** -- the inboard *fixed* mobilizer frame, rigidly attached to the parent `P`.
- **Frame_M** -- the outboard *moving* mobilizer frame, rigidly attached to the child body `B`.
- **Frame_B** -- the child body frame, anchored at the body's root (inboard) atom; every atom of the body is stored as a constant *station* `r_a` in `B`.

Two transforms are constant (built once at setup from the Z-matrix geometry) and one carries the joint motion:

- `X_PF` : `Frame_P` -> `Frame_F`, the inboard frame's placement on the parent (constant)
- `X_BM` : `Frame_B` -> `Frame_M`, the outboard frame's placement on the child body (constant)
- `X_FM(q)` : `Frame_F` -> `Frame_M`, the mobilizer transform, the ONLY configuration-dependent factor (joint DOF)

The body's pose in Ground is the ordered composition `X_GB = X_GP * X_PF * X_FM(q) * X_MB` with `X_MB = (X_BM)^-1`.

Atom positions follow as `r_a(Ground) = X_GB * r_a`. Only `X_FM(q)` changes as `q` moves; `X_PF` and `X_BM` are fixed at construction. For a `Pin` (torsion) joint the rotation is about the mobilizer Z axis, so the setup bakes a fixed bond->Z alignment (a -90 deg rotation about Y in this engine) into `X_PF`/`X_BM`, placing Z along the rotated bond.

For joints with a quaternion-parameterized rotation (`Free` and `Ball`), the four orientation coordinates are the quaternion `q = (qw, qx, qy, qz)` representing `R_FM`, and the three rotational generalized speeds are the angular velocity `w_FM = (wx, wy, wz)` of `Frame_M` relative to `Frame_F` **expressed in the parent F frame** (Simbody convention). The rotation map, the quaternion derivative, its inverse, and the second derivative are exactly Simbody's (Rotation.h: calcUnnormalizedNForQuaternion / NInv):

```text
R_FM(q) = [ 1-2(qy^2+qz^2)    2(qx*qy - qw*qz)   2(qx*qz + qw*qy)
            2(qx*qy + qw*qz)  1-2(qx^2+qz^2)     2(qy*qz - qw*qx)
            2(qx*qz - qw*qy)  2(qy*qz + qw*qx)   1-2(qx^2+qy^2) ]      (q normalized)

qdot = N(q) * w_FM ,    N(q) = (1/2) * [ -qx  -qy  -qz       (4x3)
                                          qw   qz  -qy
                                         -qz   qw   qx
                                          qy  -qx   qw ]

w_FM = NInv(q) * qdot , NInv(q) = 2 * [ -qx   qw  -qz   qy    (3x4)
                                        -qy   qz   qw  -qx
                                        -qz  -qy   qx   qw ]

qddot = N(q) * wdot_FM + N(qdot) * w_FM      (N is linear and homogeneous in its first argument)
```

Component form of the derivative -- the form the integrator and the SHAKE `q`-update must use:

```text
qdot_w = (1/2) * ( -qx*wx - qy*wy - qz*wz )
qdot_x = (1/2) * (  qw*wx + qz*wy - qy*wz )
qdot_y = (1/2) * ( -qz*wx + qw*wy + qx*wz )
qdot_z = (1/2) * (  qy*wx - qx*wy + qw*wz )
```

This is the parent-frame map. The body-frame map `qdot = (1/2) q (x) (0, w)` negates the off-diagonal (cross-coupling) terms of `N` -- giving `qdot_x = (1/2)(qw*wx - qz*wy + qy*wz)`, and likewise for `qdot_y`, `qdot_z` -- and expects `w` expressed in `Frame_M`. It **must not** be used here, because `w_FM` is expressed in `F` and `R_FM(q)` above is the standard map: pairing the standard `R_FM` with the body-frame `N` integrates the orientation with the wrong-handed angular velocity.

(Implementation status: the engine implements the parent-frame map above via a single shared helper `quaternionDotFromAngVel` used by both `Rotation::convertAngVelToQuaternionDot` and `Quat::angVelToQdot`, and the exp-map advance uses the matching left Hamilton product. An earlier revision used the body-frame map and the right product -- conflicts C-4 and C-3, the cause of the KE pump in the integrator section -- now corrected and pinned by the primitive test suite.)

## Construction of the mobilizer frames from bond geometry

`X_FM(q)` is the joint DOF (`Rot(q, Z)` for a `Pin`, identity for a `Weld`, `Translate(q)` for a Translation, quaternion+translation for a `Free`). The two **constant** frames `X_PF` and `X_BM` are rebuilt once per configuration transfer from the actual Cartesian geometry, in two stages.

### Stage 1 -- a frame per atom, from bond vectors

Each atom a gets a Ground-frame pose `frame[a]` and a cross-bond transform `xpc[a]`, built from four atoms: the atom (self), its parent, its
grandparent, and a reference child (the four atoms that define the atom's internal coordinates). With `pSelf`, `pParent`, `pGparent`, `pRefChild` their Cartesian positions:

```text
R_G : X axis = unit(pSelf - pParent)                          // along the bond, parent -> self
      Y axis = component of unit(pGparent - pParent) perp to X   (Gram-Schmidt)
      Z axis = X x Y                                          // normal to the (self,parent,gparent) plane

G        = ( R_G , pParent )                                  // a frame anchored at the PARENT atom
theta    = dihedral(pGparent, pParent, pSelf, pRefChild)      // torsion about the parent->self bond
d        = |pSelf - pParent|                                  // bond length
B        = Rot(theta, X) * Translate(d along X) * Rot(180 deg, Y)
frame[a] = G * B * Rot(180 deg, X)
xpc[a]   = B
```

So the **bond vector is the local X axis**, the **grandparent fixes the perpendicular** (the (self,parent,gparent) plane; its normal is Z), the **dihedral** rides on X, and the **bond length** is the X-translation. The two 180-degree flips (about Y inside B, then about X) reproduce the bond-center "about-face" convention of the original Molmodel build -- the parent's outboard bond-center direction and the child's inboard one point in opposite senses along the bond -- so `frame[a]` matches the original bond-center frame protocol to ~1e-14. That equivalence is exactly what the angle/chirality reconstruction cross-check verifies: rebuilding each atom's local bond-center directions from bond angles + chirality and comparing against the position-derived frame localizes
any mismatch to the offending atom. A **root** atom gets frame = (identity rotation, pSelf), xpc = identity. An atom lacking a full (parent, grandparent) lineage falls back to X = bond direction with xpc = Translate(d, X) * Rot(180 deg, Y) and no dihedral.

### Stage 2 -- mobilizer frames from the per-atom frame

For body b with root (inboard) atom `root` and parent body p with root atom `proot`:

```text
T_X_B = frame[root]                                          // the body frame B == its root atom's frame
if p is Ground:
    X_PF[b] = T_X_B                                          // Frame_F is the root-atom frame, in Ground
    X_BM[b] = identity                                       // Frame_M coincides with B
else:
    Proot_X_root = inverse(frame[proot]) * frame[root]       // pose of B's root atom in p's root-atom frame
    X_BM[b] = xpc[root] * X_to_Z                             // Frame_B -> Frame_M
    X_PF[b] = Proot_X_root * X_BM[b]                         // Frame_P -> Frame_F
where   X_to_Z = Rot(-90 deg, Y)                             // carries the bond (X) onto the mobilizer Z axis
```

`X_to_Z` is the only joint-convention-specific piece: a `Pin` rotates about the mobilizer **Z** axis, but Stage 1 puts every bond on **X**, so the -90-degree rotation about Y carries X onto Z (X -> +Z, Z -> -X). The `Frame_M` frame's Z axis therefore lies along the inboard bond, and `X_FM(q) = Rot(q, Z)` turns the body about that bond. The reference dihedral is baked into `xpc[root]` (hence into `X_BM` and `X_PF`), so at `q=0` the chain reproduces the **incoming** geometry exactly, and the `Pin` coordinate `q` is the torsion measured from that carried-over reference.

With these definitions and `X_FM(0) = Identity`, the composition must return the body frame:

```text
X_GB(0) = X_GP * X_PF * X_FM(0) * X_MB
        = frame[proot] * (Proot_X_root * X_BM) * I * (X_BM)^-1
        = frame[proot] * Proot_X_root
        = frame[proot] * inverse(frame[proot]) * frame[root]
        = frame[root] = T_X_B.
```

Any change to the frame build must preserve this identity; comparing X_GB at q = 0 against frame[root] for every body is the fastest regression test, and -- unlike an energy check -- it catches the X_MB-vs-X_BM and bond->Z sign errors.

The three factors of the chain are recomputed on deliberately different schedules, and that split is the point of the design:

- `X_PF (P->F), X_BM (B->M, hence X_MB = ~X_BM)`: once per Gibbs-block transfer, plus once at world construction. Role: the inboard bond's frozen geometry (length + angle + reference dihedral); held constant for the block's entire MD trajectory.
- `X_FM(q)`: every MD step. Role: the live joint DOF; for a `Pin`, `Rot(q, Z)` with the current torsion `q`.
- `X_GP` (the parent body's X_GB) and `X_GB` itself: every MD step. Role: the running composition; `X_GP` is constant only for a root body (parent = Ground)

So `X_PF` and `X_BM` carry the per-block **frozen** geometry while `X_FM`, `X_GP`, and `X_GB` carry the **per-step** motion. **The same transfer that rebuilds `X_PF` and X_BM `also` resets `q=0`**, making the incoming configuration the block's reference. Within a block only `X_FM(q)` and the ancestor poses move, so `X_GB(q)` is the only thing the integrator recomputes each step.

The frame graph that Stage 1 consumes is filled from the molecular bond graph with ring-closing bonds removed (parent = the bonded atom toward the root,
grandparent = its parent, refChild = the atom's first child). That acyclic tree is the same BAT / Z-matrix internal-coordinate tree the topology builder produces as `(z_i, z_j, z_k, z_l)`. The correspondence is:

- Tree topology: each atom's parent `z_j` and grandparent `z_k` is identical to the Z-matrix.
- The bond length `d = |self - parent|` and the bond angle (the grandparent fixes the frame's perpendicular) are the Z-matrix bond/angle internal coordinates, and they are frozen into `X_PF` and `X_BM`.
- Torsions are the live coordinates: a flexible dihedral about a tree bond becomes a `Pin` joint, so a Z-matrix torsion is `q = X_FM(q)`. Its reference value is baked into `X_BM` via `xpc[root]` (which carries `Rot(theta, X)`), and `X_FM(q) = Rot(q, Z)` turns about that bond.

Two qualifications are load-bearing and must not be relaxed into an identity. First, the frame kernel derives this tree and all geometry from the Cartesian positions; it shares the *tree* with the separate Z-matrix builder, not the code, and it uses the atom's first child as the local dihedral reference (pointing away from the root), whereas a textbook Z-matrix `z_l` references toward the root -- both describe the torsion about the same bond, but with a different fourth atom. Second, the frozen bond lengths and angles are not hard-constrained to idealized values: they sit at whatever the incoming geometry holds and are merely held fixed for the block (Gibbs conditioning), never a permanent constraint.

For body `b` with root atom `root = bodyRootAtom[b]` (the code's "inboard (root) atom of this body"), the inboard bond is the tree edge `(par[root], root)`:

- The **inboard** atom `A = par[root]` lives in the parent body. The `F` frame `X_PF` is fixed to the parent at this joint.
- The **outboard** atom `B = root` is the child body's root atom. The `M` frame `X_BM` and the body frame `B` are anchored here.

Which physical atom is inboard versus outboard is set per molecule by the Z-matrix root, not by any fixed convention such as N-before-CA: for a given molecule it is determined, but it must be read off that molecule's tree rather than assumed.

## The multibody dynamics engine (articulated-body recursions)

Every kinematic and mass-matrix quantity is computed by `O(n)` recursions over the n bodies of the tree via **Featherstone's articulated-body algorithm**, in spatial-operator form. The global mass matrix `M`, its inverse `M^-1`, its square root, and its determinant are **never assembled as dense matrices**. Each is applied as an operator by one or two sweeps. This linear cost is the reason the internal-coordinate representation exists, and it is what makes per-block sampling scale to large systems.

Each body `b` carries, in the Ground frame, a spatial velocity `V_GB = (w ; v)` (angular **over** linear), a spatial force `F = (torque ; force)` (same ordering), and a spatial inertia `Mk_G[b]` (the body's rigid-body spatial inertia about its own origin `Bo`, re-expressed in Ground; mass, first moment `mass*c`, and the inertia tensor). The spatial inner product is the straight component pairing `< V, F > = w * torque + v * force`, so a generalized force is `tau = H^T F` under that pairing. Four operators relate a body to its parent p:

- **Phi (rigid shift / transmission).** With `r=p_PB` the vector from the parent origin `Po` to the body origin `Bo` (in Ground), `Phi` is the 6x6 built from 3x3 blocks `[[1, r_x],[0, 1]]`, `r_x` the cross-product matrix of `r`. `~Phi` (its transpose) shifts velocities and accelerations **outward** (parent -> child); `Phi` shifts forces and inertias **inward** (child -> parent). It carries no mass - pure geometry.
- **H (hinge / joint map).** Maps a body's dof generalized speeds to its cross-mobilizer spatial velocity, `V_PB_G = H*u`. `H` is built in the `F` frame and rotated to Ground.
- **P (articulated body inertia).** The apparent spatial inertia a body exhibits as the free base of the subtree hanging below it. Unlike a rigid-body inertia, it is a full symmetric 6x6 (the apparent mass depends on the direction of push, and there is no single center of mass).
- **D, Ga.** The mobility-space (hinge) inertia `D = H^T * P * H` (dof x dof), its inverse `D^-1`, and the articulated **gain** `Ga = P * H * D^-1` assembled down the tree. The symbol `Ga` is deliberately distinct from the constraint Jacobian `G` and the kinetic energy `KE`. The inertia the parent feels across the mobilizer is the articulated shift `P^+ = P - Ga * D * Ga^T = P - P * H * D^-1 * H^T * P`.

**The recursions -- each a single sweep, O(N):**

1. **Position -- outward (root -> leaves).** Compose the frame chain `X_GB = X_GP * X_PF * X_FM(q) * X_MB`. Build `H` in Ground, the shift `Phi`, the body spatial inertia in Ground, and the per-atom stations `r_a(Ground) = X_GB * r_a`.
2. **Velocity -- outward.** Shift the parent's spatial velocity and add the mobilizer term `V_GB[b] = (~Phi[b]) * V_GB[parent] + H[b] * u[b]`.
   The same sweep accumulates the two velocity-dependent bias terms the integrator needs. The gyroscopic force (Ground frame; `I` the unit inertia about `Bo`, `c` the COM offset, `w=w_GB`): `b_gyro[b] = mass * ( w x (I w) ;  w x (w x c) )`.

   and the Coriolis bias acceleration, carried outward like the velocity (`w_GP`, `v_GP` the parent spatial velocity; `w_PB_G`, `v_PB_G` the cross-mobilizer velocity; `w_FM` the cross-mobilizer angular velocity in `F`; `r_MB_F = R_FM * r_MB` the `Mo` -> `Bo` bond vector expressed in `F`):

   ```text
   a_mob[b] = ( w_GP x w_PB_G ;
                w_GP x (v_GB - v_GP) + w_GP x v_PB_G + R_GF * ( w_FM x (w_FM x r_MB_F) ) )
   a_cor[b] = (~Phi[b]) * a_cor[parent] + a_mob[b]
   ```

   The centripetal term `R_GF (w_FM x (w_FM x r_MB_F))` is nonzero for every body whose `M` frame is offset from `Bo` (every non-root `Pin` body, offset = bond length); dropping it silently pumps energy. It is the velocity-side companion of the `X_MB` pitfall -- `r_MB = X_MB.p()`, never `X_BM.p()`.
3. **Articulated body inertia -- inward (leaves -> root).** Featherstone's backward pass:

   ```text
   P[b]   = Mk_G[b] + sum_{children c} Phi[c] * P^+[c] * (~Phi[c])
   D[b]   = H[b]^T P[b] H[b],   Ga[b] = P[b] H[b] D[b]^-1,   P^+[b] = P[b] - Ga[b] D[b] Ga[b]^T
   ```

   A `Weld` (0 dof) contributes nothing through its (empty) mobility space, so `P^+ = P`.
4. **Forward dynamics (M^-1) -- inward then outward.** Given a generalized force `f`, the acceleration `udot = M^-1 (f - f_bias)` is produced by an inward residual pass (the per-body bias `Z` and the mobility residual `eps = f - H^T * Z`) followed by an outward acceleration pass (`A_GB` and `udot`), using `D^-1` and `Ga`. The bias force `f_bias` is assembled from `b_gyro` and `a_cor` of recursion 2 (the gyroscopic force plus the Coriolis acceleration mapped through the articulated inertia): with that `f_bias` this is the equation of motion solved each Verlet substep; with `f_bias = 0` it is the pure operator `a = M^-1 * f` used for the constraint solver and the kinetic energy.
5. **Force -- inward.** The transpose of forward kinematics: per-atom Cartesian forces are reduced to a generalized force by accumulating each body's spatial force and sweeping child -> parent through `Phi^T`, reading `tau = H^T` (transmitted spatial force). This is exactly the `G^T` assembly the constraint solver uses.

**Mass-matrix operators -- all O(n), all matrix-free.** The same per-body factors `{Phi, H, D, Ga}`
yield, without ever forming M:

- **M^-1 * f:** Forward dynamics (recursion 4 with `f_bias = 0`).
- **sqrt(M^-1) z:** Outward sweep using the per-body `D^-1/2` in place of `D^-1`. It maps a white-noise vector `z` to a generalized-speed vector with covariance `M^-1`, so the equipartition draw is realized as `u = sqrt(R*T) * sqrt(M^-1) * z`, giving `u ~ N(0, R*T*M^-1)` and hence `p = M u ~ N(0, R*T*M)` at O(N), with det `M` never appearing in the proposal.
- **ln(det(M)) = sum_b ln(det(D_b)):** The factorization makes the determinant the product of the per-body hinge determinants - the Fixman tree term, read from the same `D_b = H_b^T * P_b * H_b`.
- **(1/2) u^T M u:** Kinetic energy, from one application of `M`.

**The operator factorization.** In spatial-operator form the mass matrix factors through the same local quantities"

```text
M       = (I + H * Phi * Ga) * D * (I + H * Phi * Ga)^T
M^(1/2) = (I + H * Phi * Ga) D^(1/2)
v       = D^(1/2) (I + H * Phi * Ga)^T * theta_dot       (the whitening map)
```

with `Ga` the gain operator and `v` the whitened velocity. Because each factor is local (a body and its parent only) `M`, `M^-1`, `M^(1/2)`, and
`ln det(M)` all follow from one or two sweeps over the tree, which is the O(n) claim. The generalized equipartition principle is then exact by construction: drawing the whitened `v` as unit Gaussian noise and applying `sqrt(M^-1)` gives `< u u^T > = R * T * M^-1`, i.e. equipartition of kinetic energy across the reduced coordinates despite the configuration-dependent metric.
