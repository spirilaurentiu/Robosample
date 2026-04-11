# import mdtraj as md
# import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import robosample
from openmm.app import *
from openmm.unit import *

from openmm import *

SEED = 42
TD_FIXMAN_NAME = "c4bond3_td_fixman"
TD_FIXMAN_DCD = TD_FIXMAN_NAME + "_" + str(SEED) + ".repl0.dcd"
WRITE_FREQ = 1  # Write every round
TIMESTEP_PS = 0.1  # ps
NUM_STEPS_PER_ROUND = 10  # 1 ps per round
TEMPERATURE_K = 300.0  # K
NUM_EQUILIBRATION_STEPS = 0
NUM_PRODUCTION_STEPS = 1

PRMTOP_PATH = "examples/c4bond3/c4bond3.prmtop"
RST7_PATH = "examples/c4bond3/c4bond3.rst7"


# def run_cartesian(num_steps):
#     name = "c4bond3_cartesian"
#     seed = 42
#     write_freq = 1

#     context = robosample.Context(
#         name=name,
#         seed=seed,
#         prmtop=prmtop_path,
#         inpcrd=rst7_path,
#         write_freq=write_freq,
#         testing=True,
#     )

#     sele = context.create_torsional_bonds([(2, 3)])
#     context.addTorsionalWorld(sele).addSampler(
#         timeStep=TIMESTEP, mdSteps=MD_STEPS, boostMDSteps=MD_STEPS
#     )

#     context.initialize([300.0])
#     context.RunREX(0, num_steps)


def run_td(name: str) -> None:
    context = robosample.Context(
        name=name,
        seed=SEED,
        prmtop=PRMTOP_PATH,
        inpcrd=RST7_PATH,
        write_freq=WRITE_FREQ,
        testing=True,
    )

    # Spin the middle bond
    sele = context.create_torsional_bonds([(1, 2)])
    context.addTorsionalWorld(sele).addSampler(
        timeStep=TIMESTEP_PS,
        mdSteps=NUM_STEPS_PER_ROUND,
        boostMDSteps=NUM_STEPS_PER_ROUND,
        useFixmanPotential=False,
        acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    )

    # Create one replica at the target temperature (no temperature ladder needed for testing)
    context.initialize([TEMPERATURE_K])
    context.RunREX(NUM_EQUILIBRATION_STEPS, NUM_PRODUCTION_STEPS)


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


def plot_dihedral_histogram(dcd_path, atom_indices, output_file="histogram.png"):
    traj = md.load(dcd_path, top=PRMTOP_PATH)
    indices = np.array([atom_indices])
    radians = md.compute_dihedrals(traj, indices)
    degrees = np.rad2deg(radians).flatten()

    plt.figure(figsize=(8, 6))
    plt.hist(
        degrees,
        bins=50,  # 7.2 degree increments.
        range=[-180, 180],
        density=True,  # Plot PDF instead of counts
        color="skyblue",
        edgecolor="black",
        alpha=0.7,
    )

    plt.xlabel(r"$\alpha$ (degrees)")
    plt.ylabel(r"$\rho(\alpha)$")
    plt.xticks([-180, -90, 0, 90, 180])
    plt.axhline(y=1 / 360, color="black", linestyle="--", label=r"$\frac{1}{2\pi}$")

    plt.tight_layout()
    # plt.savefig(output_file, dpi=300)
    plt.show()
    plt.close()


def test_fixman_potential():
    run_td(TD_FIXMAN_NAME)
    # plot_dihedral_histogram(
    #     TD_FIXMAN_DCD,
    #     (0, 1, 2, 3),
    #     output_file="c4bond3_td_fixman_histogram.png",
    # )
