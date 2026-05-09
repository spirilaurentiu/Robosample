import mdtraj as md
import numpy as np
import pytest
import robosample


# Robosample will put the original conformation in the frame of the trajectory
# This happens before integration starts, so we are testing topology and world creation, coordinate transfer and DCD writing
@pytest.fixture(scope="module")
def simulation_results():
    name = "1APQ.test.rigid"
    seed = 42
    prmtop = "examples/1APQ.prmtop"
    xyz = "examples/1APQ.rst7"
    write_freq = 0
    equil_steps = 0
    prod_steps = 0

    context = robosample.Context(
        name=name,
        seed=seed,
        prmtop=prmtop,
        inpcrd=xyz,
        write_freq=write_freq,
        testing=True,
    )

    # Put torsions on middle bonds in standard dihedrals and integrate with 10 fs time step for 1 ps
    # Torsions are protein phi, psi and chi1
    sele = context.getDefaultBonds("standard")
    context.add_torsional_world(sele).add_sampler(
        timeStep=0.01, mdSteps=0, boostMDSteps=0
    )

    # Add one replica at 300 K
    context.initialize([300.0])

    # Run the simulation for no equilibration steps and 100 production steps totaling 100 ps of simulation time
    context.RunREX(equil_steps, prod_steps)

    # Load the trajectory
    traj = md.load("1APQ.test.rigid_42.repl0.dcd", top=prmtop)

    # Load the reference structure for comparison
    reference = md.load(xyz, top=prmtop)

    return context, traj, reference


def _assert_bond_constancy(simulation_results, bond_indices, tolerance=1e-5):
    context, traj, reference = simulation_results
    if len(bond_indices) == 0:
        return

    # Compute distances: (frames, nbonds)
    distances_ref = md.compute_distances(reference, bond_indices, periodic=False)[0]
    distances_traj = md.compute_distances(traj, bond_indices, periodic=False)

    # Calculate absolute deviations
    diffs = np.abs(distances_traj - distances_ref)
    max_deviations = np.max(diffs, axis=0)

    violations = []
    for i, dev in enumerate(max_deviations):
        if dev > tolerance:
            # Find the actual max value reached in the trajectory for this bond
            max_val = distances_traj[np.argmax(diffs[:, i]), i]

            a1, a2 = bond_indices[i]
            atom_names = f"{context.getAtomNameByPrmtopIndex(a1)}-{context.getAtomNameByPrmtopIndex(a2)}"
            violations.append(
                f"  - {atom_names} ({a1},{a2}): Ref={distances_ref[i]:.5f}nm, "
                f"MaxVal={max_val:.5f}nm, MaxDev={dev:.2e}nm, ring closing={context.getWorld(0).isBondRingClosing(a1, a2)}"
            )

    if violations:
        header = f"Bond constraints violated (Tol: {tolerance}nm):"
        pytest.fail(f"{header}\n" + "\n".join(violations))


def _assert_angle_constancy(simulation_results, angle_indices, tolerance=1e-5):
    context, traj, reference = simulation_results
    if len(angle_indices) == 0:
        return

    # Compute angles: (frames, nangles)
    angles_ref = md.compute_angles(reference, angle_indices, periodic=False)[0]
    angles_traj = md.compute_angles(traj, angle_indices, periodic=False)

    # Calculate absolute deviations
    diffs = np.abs(angles_traj - angles_ref)
    max_deviations = np.max(diffs, axis=0)

    violations = []
    for i, dev in enumerate(max_deviations):
        if dev > tolerance:
            # Find the actual max value reached in the trajectory for this angle
            max_val = angles_traj[np.argmax(diffs[:, i]), i]

            a1, a2, a3 = angle_indices[i]
            atom_names = f"{context.getAtomNameByPrmtopIndex(a1)}-{context.getAtomNameByPrmtopIndex(a2)}-{context.getAtomNameByPrmtopIndex(a3)}"
            violations.append(
                f"  - {atom_names} ({a1},{a2},{a3}): Ref={angles_ref[i]:.5f}rad, "
                f"MaxVal={max_val:.5f}rad, MaxDev={dev:.2e}rad"
            )

    if violations:
        header = f"Angle constraints violated (Tol: {tolerance}rad):"
        pytest.fail(f"{header}\n" + "\n".join(violations))


def _assert_torsion_constancy(simulation_results, torsion_indices, tolerance=1e-5):
    context, traj, reference = simulation_results
    if len(torsion_indices) == 0:
        return

    tors_ref = md.compute_dihedrals(reference, torsion_indices, periodic=False)[0]
    tors_traj = md.compute_dihedrals(traj, torsion_indices, periodic=False)

    # Calculate circular differences: ensures -pi and pi are treated as the same point
    raw_diffs = tors_traj - tors_ref
    circular_diffs = np.abs(np.arctan2(np.sin(raw_diffs), np.cos(raw_diffs)))

    # max_deviations: (nangles,) -> the furthest each torsion moved from ref
    max_deviations = np.max(circular_diffs, axis=0)

    violations = []

    for i, dev in enumerate(max_deviations):
        a1, a2, a3, a4 = torsion_indices[i]
        atom_names = f"{context.getAtomNameByPrmtopIndex(a1)}-{context.getAtomNameByPrmtopIndex(a2)}-{context.getAtomNameByPrmtopIndex(a3)}-{context.getAtomNameByPrmtopIndex(a4)}"

        if dev > tolerance:
            max_val = tors_traj[np.argmax(circular_diffs[:, i]), i]
            violations.append(
                f"  - {atom_names} ({a1},{a2},{a3},{a4}): Ref={tors_ref[i]:.5f}, "
                f"MaxVal={max_val:.5f}, MaxDev={dev:.2e} > Tol:{tolerance}"
            )

    if violations:
        header = "Torsion behavior check failed:"
        pytest.fail(f"{header}\n" + "\n".join(violations))


def test_num_frames(simulation_results):
    context, traj, reference = simulation_results
    assert traj.n_frames == 1, f"Expected 1 frames, but got {traj.n_frames}"


def test_transfer_coordinates(simulation_results):
    context, traj, reference = simulation_results
    for error in context.getWorld(0).get_coordinate_transfer_errors():
        for residual in error.matchResiduals:
            assert residual < 1e-5, (
                f"Coordinate transfer residual too high: {residual:.2e} nm"
            )
        assert error.cartesian < 1e-5, (
            f"Cartesian coordinate transfer error too high: {error.cartesian:.2e} nm"
        )
        assert error.cartesianMax < 1e-5, (
            f"Max Cartesian coordinate transfer error too high: {error.cartesianMax:.2e} nm"
        )
        assert error.bonds < 1e-5, (
            f"Bond length transfer error too high: {error.bonds:.2e} nm"
        )
        assert error.bondsMax < 1e-5, (
            f"Max bond length transfer error too high: {error.bondsMax:.2e} nm"
        )
        assert error.angles < 1e-5, (
            f"Angle transfer error too high: {error.angles:.2e} rad"
        )
        assert error.anglesMax < 1e-5, (
            f"Max angle transfer error too high: {error.anglesMax:.2e} rad"
        )
        assert error.properDihedrals < 1e-5, (
            f"Proper dihedral transfer error too high: {error.properDihedrals:.2e} rad"
        )
        assert error.properDihedralsMax < 1e-5, (
            f"Max proper dihedral transfer error too high: {error.properDihedralsMax:.2e} rad"
        )
        assert error.improperDihedrals < 1e-5, (
            f"Improper dihedral transfer error too high: {error.improperDihedrals:.2e} rad"
        )
        assert error.improperDihedralsMax < 1e-5, (
            f"Max improper dihedral transfer error too high: {error.improperDihedralsMax:.2e} rad"
        )


def test_rigid_bonds(simulation_results):
    context, traj, reference = simulation_results
    bonds = [
        [bond.prmtop_indices[0], bond.prmtop_indices[1]]
        for bond in context.bond_stretches
    ]
    _assert_bond_constancy(simulation_results, bonds)


def test_rigid_angles(simulation_results):
    context, traj, reference = simulation_results
    angles = [
        [angle.prmtop_indices[0], angle.prmtop_indices[1], angle.prmtop_indices[2]]
        for angle in context.bond_bends
    ]
    _assert_angle_constancy(simulation_results, angles)


def test_torsions(simulation_results):
    context, traj, reference = simulation_results
    periodic_torsions = [
        [
            torsion.prmtop_indices[0],
            torsion.prmtop_indices[1],
            torsion.prmtop_indices[2],
            torsion.prmtop_indices[3],
        ]
        for torsion in context.periodic_torsions
    ]
    improper_harmonic_torsions = [
        [
            torsion.prmtop_indices[0],
            torsion.prmtop_indices[1],
            torsion.prmtop_indices[2],
            torsion.prmtop_indices[3],
        ]
        for torsion in context.improper_harmonic_torsions
    ]
    torsions = periodic_torsions + improper_harmonic_torsions
    _assert_torsion_constancy(simulation_results, torsions)
