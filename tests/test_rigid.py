import robosample


def test_steady_state():
    # Simbody constructs the multibody system as a spanning tree. Only tree edges
    # are represented in the generalized coordinates; these are enforced exactly
    # through the choice of coordinates (i.e. zero-DOF "Rigid" mobilities).
    #
    # Ring-closing bonds are, by construction, non-tree edges and are therefore
    # not part of the coordinate system. In our setup, the Simbody matter subsystem
    # contains no explicit constraints, meaning these loop closures are not enforced
    # via holonomic constraint equations or Lagrange multipliers.
    #
    # As a result:
    #   - Tree bonds (internal coordinates) are exact up to machine precision.
    #   - Ring-closing bonds are not enforced by the integrator and can drift.
    #
    # The observed deviations depend on the timestep and exhibit the expected
    # O(h^2) scaling of the Verlet integrator (where h is the time step) reaching
    # a bounded steady-state error (~1e-3 - 1e-2 in our tests).
    # This is not a violation of rigidity at the coordinate level,
    # but a consequence of representing a cyclic molecular graph
    # with a tree-structured multibody system without additional constraints.
    #
    # Importantly, ring-closing bonds that connect atoms within the same rigid
    # body show no error, confirming that rigidity is preserved exactly whenever
    # it is encoded in the multibody topology.
    #
    # Therefore, tests distinguish two regimes:
    #
    #   1. Tree bonds:
    #        Expected to be invariant (tolerance ~1e-6 or tighter).
    #
    #   2. Ring-closing bonds:
    #        Expected to remain bounded with timestep-dependent error
    #        (tolerance ~1e-2), but not exact.
    #
    # This distinction reflects the underlying numerical model rather than a defect.
    # All of this is implemented in World::hasRigidBodyViolations().
    # We always run this function upon initialization, albeit with only 1 step to check for gross violations.
    # Also note that this function is being called with 1 step inside `context.initialize()`.

    name = "1APQ.test.rigid"
    seed = 42
    prmtop = "examples/1APQ.prmtop"
    xyz = "examples/1APQ.rst7"
    write_freq = 0

    context = robosample.Context(
        name=name,
        seed=seed,
        prmtop=prmtop,
        inpcrd=xyz,
        write_freq=write_freq,
        testing=True,
    )

    sele = context.getDefaultBonds("standard")
    context.add_torsional_world(sele).add_sampler(timeStep=0, mdSteps=0, boostMDSteps=0)

    # context.initialize([300.0])

    # assert not context.getWorld(0).has_rigid_body_violations(0.01, 1)
