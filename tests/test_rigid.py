import robosample
import numpy as np
import mdtraj as md
import pytest

name = 'ala-dipeptide.test.rigid'
seed = 42
prmtop = 'examples/ala-dipeptide.prmtop'
inpcrd = 'examples/ala-dipeptide.rst7'
write_freq = 1
equil_steps = 0
prod_steps = 10

context = robosample.Context(name=name, seed=seed, prmtop=prmtop, inpcrd=inpcrd, write_freq=write_freq, testing=True)

# Put torsions on middle bonds in standard dihedrals and integrate with 10 fs time step for 1 ps
sele = context.getDefaultBonds('standard')
context.addTorsionalWorld(sele).addSampler(timeStep=0.01, mdSteps=100, boostMDSteps=100)

# Add one replica at 300 K
context.initialize([300.0])

# Run the simulation for no equilibration steps and 100 production steps totaling 100 ps of simulation time
context.RunREX(equil_steps, prod_steps)

# Load the trajectory
reference = md.load(inpcrd, top=prmtop)
traj = md.load('ala-dipeptide.test.rigid_42.repl0.dcd', top=prmtop)

def _assert_bond_constancy(bond_indices, tolerance=1e-4, label="Bond"):
    if len(bond_indices) == 0:
        return

    # Compute distances: (frames, nbonds)
    distances_ref = md.compute_distances(reference, bond_indices)[0]
    distances_traj = md.compute_distances(traj, bond_indices)
    
    # Vectorized max deviation calculation
    max_deviations = np.max(np.abs(distances_traj - distances_ref), axis=0)
    
    violations = []
    for i, dev in enumerate(max_deviations):
        if dev > tolerance:
            a1, a2 = bond_indices[i]
            atom_names = f"{context.getAtomName(a1)}-{context.getAtomName(a2)}"
            violations.append(
                f"  - {atom_names} ({a1},{a2}): Ref={distances_ref[i]:.5f}nm, MaxDev={dev:.2e}nm"
            )

    if violations:
        header = f"{label} constraints violated (Tol: {tolerance}nm):"
        pytest.fail(f"{header}\n" + "\n".join(violations))

def _assert_angle_constancy(angle_indices, tolerance=1e-4, label="Angle"):
    if len(angle_indices) == 0:
        return

    # Compute angles: (frames, nangles)
    angles_ref = md.compute_angles(reference, angle_indices)[0]
    angles_traj = md.compute_angles(traj, angle_indices)
    
    # Vectorized max deviation calculation
    max_deviations = np.max(np.abs(angles_traj - angles_ref), axis=0)
    
    violations = []
    for i, dev in enumerate(max_deviations):
        if dev > tolerance:
            a1, a2, a3 = angle_indices[i]
            atom_names = f"{context.getAtomName(a1)}-{context.getAtomName(a2)}-{context.getAtomName(a3)}"
            violations.append(
                f"  - {atom_names} ({a1},{a2},{a3}): Ref={angles_ref[i]:.5f}rad, MaxDev={dev:.2e}rad"
            )

    if violations:
        header = f"{label} constraints violated (Tol: {tolerance}rad):"
        pytest.fail(f"{header}\n" + "\n".join(violations))

def test_rigid_bonds():
    bonds = context.getWorld(0).getTestRigidBonds()
    _assert_bond_constancy(bonds, label="Rigid Bond")

def test_non_rigid_bonds():
    bonds = context.getWorld(0).getTestNonRigidBonds()
    _assert_bond_constancy(bonds, label="Non-Rigid Bond")

def test_rigid_angles():
    angles = context.getWorld(0).getTestRigidAngles()
    _assert_angle_constancy(angles, label="Rigid Angle")

def test_non_rigid_angles():
    angles = context.getWorld(0).getTestNonRigidAngles()
    _assert_angle_constancy(angles, label="Non-Rigid Angle")