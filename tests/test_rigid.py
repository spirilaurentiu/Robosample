import robosample
import numpy as np
import mdtraj as md
import pytest

name = '1APQ.test.rigid'
seed = 42
prmtop = 'examples/1APQ.prmtop'
xyz = 'examples/1APQ.rst7'
write_freq = 1
equil_steps = 0
prod_steps = 100

context = robosample.Context(name=name, seed=seed, prmtop=prmtop, inpcrd=xyz, write_freq=write_freq, testing=True)

# Put torsions on middle bonds in standard dihedrals and integrate with 10 fs time step for 1 ps
# Torsions are protein phi, psi and chi1
sele = context.getDefaultBonds('standard')
context.addTorsionalWorld(sele).addSampler(timeStep=0.01, mdSteps=100, boostMDSteps=100)

# Add one replica at 300 K
context.initialize([300.0])

# Run the simulation for no equilibration steps and 100 production steps totaling 100 ps of simulation time
context.RunREX(equil_steps, prod_steps)

# Load the trajectory
traj = md.load('1APQ.test.rigid_42.repl0.dcd', top=prmtop)
# reference = md.load(xyz, top=prmtop)
reference = traj[0]

def _assert_bond_constancy(bond_indices, tolerance=1e-5, label="Bond"):
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

            import matplotlib.pyplot as plt
            plt.plot(distances_traj[:, i], label=f"Bond {i}")
            plt.axhline(distances_ref[i], color='red', linestyle='--', label='Reference')
            plt.title(f"Bond {i} distance over time")
            plt.xlabel("Frame")
            plt.ylabel("Distance (nm)")
            plt.legend()
            plt.show()
            
            a1, a2 = bond_indices[i]
            atom_names = f"{context.getAtomNameByPrmtopIndex(a1)}-{context.getAtomNameByPrmtopIndex(a2)}"
            violations.append(
                f"  - {atom_names} ({a1},{a2}): Ref={distances_ref[i]:.5f}nm, "
                f"MaxVal={max_val:.5f}nm, MaxDev={dev:.2e}nm, ring closing={context.getWorld(0).isBondRingClosing(a1, a2)}"
            )

    if violations:
        header = f"{label} constraints violated (Tol: {tolerance}nm):"
        pytest.fail(f"{header}\n" + "\n".join(violations))

def _assert_angle_constancy(angle_indices, tolerance=1e-5, label="Angle"):
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
        header = f"{label} constraints violated (Tol: {tolerance}rad):"
        pytest.fail(f"{header}\n" + "\n".join(violations))

def _assert_torsion_constancy(torsion_indices, mode, tolerance=1e-5, label="Torsion"):
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
        
        if mode == "constant":
            if dev > tolerance:
                max_val = tors_traj[np.argmax(circular_diffs[:, i]), i]
                violations.append(
                    f"  - {atom_names} ({a1},{a2},{a3},{a4}): Ref={tors_ref[i]:.5f}, "
                    f"MaxVal={max_val:.5f}, MaxDev={dev:.2e} > Tol:{tolerance}"
                )
        
        elif mode == "dynamic":
            if dev < tolerance:
                violations.append(
                    f"  - {atom_names} ({a1},{a2},{a3},{a4}): Static! "
                    f"Max deviation was only {dev:.5f}rad (Min required: {tolerance})"
                )
        else:
            raise ValueError("mode must be 'constant' or 'dynamic'")

    if violations:
        header = f"{label} behavior [{mode}] check failed:"
        pytest.fail(f"{header}\n" + "\n".join(violations))

def test_rigid_bonds():
    # In torsional dynamics, all bond lengths should be constant, including non-rigid ones since they are not allowed to stretch
    rigid_bonds = context.getWorld(0).getTestRigidBonds()
    non_rigid_bonds = context.getWorld(0).getTestNonRigidBonds()

    # These indices are returned as pairs of prmtop indices
    all_bonds = rigid_bonds + non_rigid_bonds
    ref_bonds = [(b[0].index, b[1].index) for b in traj.topology.bonds]

    normalized_all = {tuple(sorted(bond)) for bond in all_bonds}
    normalized_ref = {tuple(sorted(bond)) for bond in ref_bonds}
    
    assert normalized_all == normalized_ref, f"Bonds do not match! Missing: {normalized_ref - normalized_all} Extra: {normalized_all - normalized_ref}"

    _assert_bond_constancy(rigid_bonds, label="Rigid Bond")
    _assert_bond_constancy(non_rigid_bonds, label="Non-Rigid Bond")

# def test_rigid_angles():
#     rigid_angles = context.getWorld(0).getTestRigidAngles()
#     non_rigid_angles = context.getWorld(0).getTestNonRigidAngles()

#     _assert_angle_constancy(rigid_angles, label="Rigid Angle")
#     _assert_angle_constancy(non_rigid_angles, label="Non-Rigid Angle")

# def test_rigid_proper_torsions():
#     rigid_torsions = context.getWorld(0).getTestRigidProperTorsions()
#     non_rigid_torsions = context.getWorld(0).getTestNonRigidProperTorsions()

#     _assert_torsion_constancy(rigid_torsions, mode="constant", label="Rigid Proper Torsion")
#     _assert_torsion_constancy(non_rigid_torsions, mode="dynamic", label="Non-Rigid Proper Torsion")

# def test_rigid_improper_torsions():
#     rigid_torsions = context.getWorld(0).getTestRigidImproperTorsions()
#     non_rigid_torsions = context.getWorld(0).getTestNonRigidImproperTorsions()
    
#     _assert_torsion_constancy(rigid_torsions, mode="constant", label="Rigid Improper Torsion")
#     _assert_torsion_constancy(non_rigid_torsions, mode="dynamic", label="Non-Rigid Improper Torsion")
