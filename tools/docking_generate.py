import argparse
import mdtraj as md
import numpy as np
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import AllChem
import subprocess
import tempfile
import os
from scipy.spatial.transform import Rotation as R
import concurrent.futures

import openmm as mm
import openmm.app as app
import openmm.unit as unit
from sys import stdout
import argparse
import parmed as pmd

import flexor
import mdtraj as md
import argparse
import robosample

def point_in_hull(point, hull, tolerance=1e-12):
    return all((np.dot(eq[:-1], point) + eq[-1] <= tolerance) for eq in hull.equations)

def sample_points_in_convex_hull_volume(points, num_samples):
    # Create the original convex hull
    original_hull = ConvexHull(points)

    # Calculate normal vectors for each face
    normal_vectors = original_hull.equations[:, :-1]
    
    # Normalize the normal vectors
    normal_vectors = normal_vectors / np.linalg.norm(normal_vectors, axis=1)[:, np.newaxis]
    
    # Create the inner convex hull
    inner_distance = 1
    inner_points = np.array([p + inner_distance * normal_vectors for p in points]).reshape(-1, 3)
    inner_hull = ConvexHull(inner_points)

    # Create the outer convex hull
    outer_distance = 2
    outer_points = np.array([p + outer_distance * normal_vectors for p in points]).reshape(-1, 3)
    outer_hull = ConvexHull(outer_points)
    
    # Initialize an empty list to store valid samples
    valid_samples = []
    
    while len(valid_samples) < num_samples:
        # Generate a random point within the bounding box of the outer hull
        random_point = np.random.uniform(outer_hull.min_bound, outer_hull.max_bound)
        
        # Check if the point is inside the outer hull but outside the inner hull
        if point_in_hull(random_point, outer_hull) and not point_in_hull(random_point, inner_hull):
            valid_samples.append(random_point)
    
    return np.array(valid_samples)

def translate_mol2(i, input_mol2_file, point, output_mol2_file):
    # Load the mol2 file
    mol = Chem.MolFromMol2File(input_mol2_file, removeHs=False)
    
    # Get the conformer
    conf = mol.GetConformer()

    # Periodic table to get atomic masses
    periodic_table = Chem.GetPeriodicTable()

    # Calculate the center of mass (weighted by atomic masses)
    center_of_mass = np.zeros(3)
    total_mass = 0.0
    num_atoms = mol.GetNumAtoms()
    
    for atom_idx in range(num_atoms):
        atom = mol.GetAtomWithIdx(atom_idx)
        pos = conf.GetAtomPosition(atom_idx)
        
        # Get the mass of the atom using its atomic number
        mass = periodic_table.GetAtomicWeight(atom.GetAtomicNum())
        center_of_mass += mass * np.array([pos.x, pos.y, pos.z])
        total_mass += mass
    
    center_of_mass /= total_mass

    # Center the atoms to the origin using the center of mass
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        centered_pos = np.array([pos.x, pos.y, pos.z]) - center_of_mass
        conf.SetAtomPosition(atom_idx, centered_pos)

    # Generate a random rotation matrix
    np.random.seed(int.from_bytes(os.urandom(4), byteorder='little'))
    random_rotation = R.random()
    rotation_matrix = random_rotation.as_matrix()

    # Apply rotation and translation to all atom positions
    for atom_idx in range(num_atoms):
        pos = conf.GetAtomPosition(atom_idx)
        centered_pos = np.array([pos.x, pos.y, pos.z])
        rotated_pos = rotation_matrix @ centered_pos
        new_pos = rotated_pos + point
        conf.SetAtomPosition(atom_idx, new_pos)
    
    # Save the modified molecule to a new mol file
    temp_mol_file = f'temp_{i}.mol'
    Chem.MolToMolFile(mol, temp_mol_file)

    # Convert the mol file to mol2
    subprocess.run(['obabel', '-imol', temp_mol_file, '-O', output_mol2_file])
    os.remove(temp_mol_file)

# sample points -> translate mol2 -> translate_i.mol2 (obabel)
# translate_0.mol2 -> amber.mol2 (antechamber) -> amber.prepi (antechamber) -> amber.frcmod (parmchk2)
# tleap

def generate_system(i, point, mol2_file, pdb_file):
    translated_mol2 = f'translated_{i}.mol2'
    amber_mol2 = f'amber_{i}.mol2'
    prepi_file = f'fps_{i}.prepi'
    frcmod_file = f'fps_{i}.frcmod'
    tleap_file_path = f'tleap_{i}.in'
    prmtop = f'rage.fps.{i}.prmtop'
    rst7 = f'rage.fps.{i}.rst7'
    minRst7 = f'rage.fps.{i}.min.rst7'

    # Convert point around protein from nm to angstrom and translate the ligand
    translate_mol2(i, mol2_file, point * 10, translated_mol2)

    # Prepare amber input
    # In multithreaded environments, we must not delete temporary files (-pf y)
    subprocess.run(['antechamber', '-i', translated_mol2, '-fi', 'mol2', '-o', amber_mol2, '-fo', 'mol2', '-c', 'gas', '-at', 'amber', '-pf', 'y'])
    subprocess.run(['antechamber', '-i', amber_mol2, '-fi', 'mol2', '-o', prepi_file, '-fo', 'prepi', '-c', 'gas', '-at', 'amber', '-pf', 'y'])
    subprocess.run(['parmchk2', '-i', amber_mol2, '-f', 'mol2', '-o', frcmod_file])

    tleap_input = f"""
    source leaprc.protein.ff19SB
    source leaprc.gaff
    loadamberprep {prepi_file}
    loadamberparams {frcmod_file}

    protein = loadpdb {pdb_file}
    ligand = loadmol2 {amber_mol2}
    complex = combine {{ligand protein}}

    saveamberparm complex {prmtop} {rst7}
    quit
    """
    # Write tleap_input to a temporary file
    with open(tleap_file_path, 'w') as f:
        f.write(tleap_input)
    subprocess.run(['tleap', '-f', tleap_file_path])

    # Minimize
    ommPrmtop = pmd.load_file(prmtop)
    ommRst7 = pmd.amber.Rst7(rst7)

    # Create an OpenMM integrator and simulation context
    system = ommPrmtop.createSystem(nonbondedMethod=app.NoCutoff, nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds, 0.002*unit.picoseconds)
    platform = mm.Platform.getPlatformByName('CUDA')  # Use 'CUDA' or 'CPU' depending on your setup
    simulation = app.Simulation(ommPrmtop.topology, system, integrator, platform)
    simulation.context.setPositions(ommRst7.positions)
    simulation.minimizeEnergy() # maxIterations=500
    minimized_positions = simulation.context.getState(getPositions=True).getPositions()
    ommPrmtop.coordinates = minimized_positions
    ommPrmtop.save(minRst7, format='rst7', overwrite=True)

    # Delete temporary files
    os.remove(translated_mol2)
    os.remove(amber_mol2)
    os.remove(prepi_file)
    os.remove(frcmod_file)
    os.remove(tleap_file_path)
    os.remove('leap.log')
    os.remove(rst7)

    return prmtop, minRst7

def simulate(i, prmtop, rst7):
    # prepare flexor generator
    mdtrajObj = md.load(rst7, top=prmtop)
    flexorObj = flexor.Flexor(mdtrajObj)

    # create robosample context
    c = robosample.Context(f"rage.fps.{i}", 42, 1, 1, robosample.RunType.REMC, 1, 0)
    c.setPdbRestartFreq(0) # WRITE_PDBS
    c.setPrintFreq(1) # PRINT_FREQ
    c.loadAmberSystem(prmtop, rst7)

    # samplers
    sampler = robosample.SamplerName.HMC # rename to type
    thermostat = robosample.ThermostatName.ANDERSEN

    # openmm cartesian
    flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
    c.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # sidechains pins
    flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
    c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # ramachandran pins
    flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
    c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # rigid all
    flex = [robosample.BondFlexibility(-1, 3466, robosample.BondMobility.Free)]
    c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    c.getWorld(0).addSampler(sampler, robosample.IntegratorName.OMMVV, thermostat, False)
    c.getWorld(1).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)
    c.getWorld(2).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)
    c.getWorld(3).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)

    accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
    timesteps = [0.001, 0.004, 0.008, 0.1]
    world_indexes = [0, 1, 2, 3]
    mdsteps = boost_md_steps = [50, 20, 20, 10]
    distort_options = [0, 0, 0, 0]
    distort_args = ["0", "0" , "0", "0"]
    flow = [0, 0, 0, 0]
    work = [0, 0, 0, 0]

    # Calculate the common ratio for geometric progression
    t_min = 300
    t_max = 600
    num_replicas = 1

    if num_replicas > 1:
        r = (t_max / t_min) ** (1 / (num_replicas - 1))
    else:
        r = 1  # When there is only one replica, r is set to 1

    # Generate the temperatures
    temperatures = [t_min * r ** i for i in range(num_replicas)]
    boost_temperatures = [t_min * r ** i for i in range(num_replicas)]  # used for openmm velocities

    for i in range(num_replicas):
        c.addReplica(i)
        c.addThermodynamicState(i, temperatures[i], boost_temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, world_indexes, timesteps, mdsteps, boost_md_steps)

    # initialize the simulation
    c.initialize()

    # start the simulation
    c.RunREXNew(10, 100)

if __name__ == "__main__":
    # parser = argparse.ArgumentParser(description="Process PDB and MOL2 files with a specified distance.")
    # parser.add_argument('pdb_file', type=str, help='Path to the PDB file (protein)')
    # parser.add_argument('mol2_file', type=str, help='Path to the MOL2 file (ligand)')

    # args = parser.parse_args()

    # pdb_file = args.pdb_file

    # echo "C1CCC(CC1)N(CC2=CC=CC=C2)C(=O)C3=CC=C(C=C3)Cl" | obabel -ismi -omol2 --gen3d > fps.mol2
    mol2_file = 'fps.mol2'
    pdb_file = 'rage.capped.noh.pdb'

    # mamba install nvidia/label/cuda-12.3.2::cuda-toolkit

    traj = md.load_pdb(pdb_file)

    pts = traj.xyz[0]
    hull = ConvexHull(pts)

    # Sample random points on the surface of the convex hull and project them outward
    num_points = 8 # Number of points to sample
    sampled_points = sample_points_in_convex_hull_volume(pts, num_points)

    shouldPlot = False
    if shouldPlot:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")

        # Plot defining corner points
        ax.plot(pts.T[0], pts.T[1], pts.T[2], "ko")

        # Plot sampled points
        ax.plot(sampled_points.T[0], sampled_points.T[1], sampled_points.T[2], "bo")

        # Plot the convex hull
        for s in hull.simplices:
            s = np.append(s, s[0])  # Here we cycle back to the first coordinate
            ax.plot(pts[s, 0], pts[s, 1], pts[s, 2], "r-")

        # Make axis label
        for i in ["x", "y", "z"]:
            eval("ax.set_{:s}label('{:s}')".format(i, i))

        plt.show()

    prmtops, rst7s = [], []
    for i, point in enumerate(sampled_points):
        p, r = generate_system(i, point, mol2_file, pdb_file)
        prmtops.append(p)
        rst7s.append(r)

    # # Use ThreadPoolExecutor to run process_point in multiple threads
    # # max_workers = os.cpu_count() # will return the number of logical cores
    # max_workers = 4
    # with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
    #     futures = [executor.submit(simulate, i, p, r) for i, (p, r) in enumerate(zip(prmtops, rst7s))]

    # # Optionally, wait for all futures to complete
    # for future in futures:
    #     try:
    #         result = future.result()
    #     except Exception as e:
    #         print(f"An exception occurred: {e}")
