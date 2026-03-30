# import robosample
# import numpy as np
# import mdtraj as md

# from openmm.app import *
# from openmm import *
# from openmm.unit import *

# prmtop_path = 'examples/c4bond3/c4bond3.prmtop'
# rst7_path = 'examples/c4bond3/c4bond3.rst7'

# def run_cartesian(num_steps):
#     name = 'c4bond3_cartesian'
#     seed = 42
#     write_freq = 1

#     context = robosample.Context(name=name, seed=seed, prmtop=prmtop_path, inpcrd=rst7_path, write_freq=write_freq, testing=True)

#     sele = context.create_torsional_bonds([(2,3)])
#     context.addCartesianWorld().addSampler(timeStep=0.001, mdSteps=10, boostMDSteps=0)

#     context.initialize([300.0])
#     context.RunREX(0, num_steps)
    
# def run_td(num_steps):
#     name = 'c4bond3_td'
#     seed = 42
#     write_freq = 1

#     context = robosample.Context(name=name, seed=seed, prmtop=prmtop_path, inpcrd=rst7_path, write_freq=write_freq, testing=True)

#     sele = context.create_torsional_bonds([(2,3)])
#     context.addTorsionalWorld(sele).addSampler(timeStep=0.001, mdSteps=10, boostMDSteps=0)

#     context.initialize([300.0])
#     context.RunREX(0, num_steps)

# def get_dihedral_values(dcd_path, prmtop_path, atom_indices):
#     """
#     Loads a DCD trajectory and computes dihedrals for a specific set of atoms.
    
#     Parameters:
#     - dcd_path: str, path to the .dcd file
#     - prmtop_path: str, path to the .prmtop (topology) file
#     - atom_indices: tuple or list of 4 integers (0-indexed)
    
#     Returns:
#     - numpy array of dihedral angles in degrees
#     """
#     # 1. Load the trajectory
#     # MDTraj needs the topology file to understand the DCD structure
#     traj = md.load(dcd_path, top=prmtop_path)
    
#     # 2. Reshape indices for MDTraj 
#     # compute_dihedrals expects a 2D array of shape (n_dihedrals, 4)
#     indices = np.array([atom_indices])
    
#     # 3. Compute dihedrals (returns values in radians)
#     radians = md.compute_dihedrals(traj, indices)
    
#     # 4. Convert to degrees and flatten to a 1D array
#     degrees = np.rad2deg(radians).flatten()
    
#     return degrees


# import mdtraj as md
# import numpy as np
# import matplotlib.pyplot as plt

# def plot_dihedral_histogram(dcd_path, prmtop_path, atom_indices, output_file='histogram.png'):
#     # 1. Load trajectory and compute dihedrals
#     traj = md.load(dcd_path, top=prmtop_path)
#     indices = np.array([atom_indices])
#     radians = md.compute_dihedrals(traj, indices)
#     degrees = np.rad2deg(radians).flatten()

#     # 2. Create the plot
#     plt.figure(figsize=(8, 6))
#     plt.hist(degrees, bins=90, range=[-180, 180], color='skyblue', edgecolor='black', alpha=0.7)
    
#     # 3. Formatting with LaTeX
#     plt.title(f'Distribution of Dihedral Angle for Atoms {atom_indices}')
#     plt.xlabel(r'Dihedral Angle $\phi$ (degrees)')
#     plt.ylabel('Frequency')
#     plt.grid(axis='y', linestyle='--', alpha=0.7)
    
#     # Set x-ticks to common angles
#     plt.xticks(np.arange(-180, 181, 60))
#     plt.xlim(-180, 180)

#     # 4. Save the plot
#     plt.tight_layout()
#     # plt.savefig(output_file, dpi=300)
#     plt.show()
#     plt.close()

# # Example usage:
# # run_openmm(1000)
# plot_dihedral_histogram('c4bond3_42.repl0.dcd', prmtop_path, (0, 1, 2, 3))