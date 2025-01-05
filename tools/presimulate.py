import argparse
from openmm import app
from openmm import *
from openmm.unit import *
from sys import stdout

def run_simulation(prmtop_file, rst7_file, dcd_file):
    # Load the topology and positions
    prmtop = app.AmberPrmtopFile(prmtop_file)
    inpcrd = app.AmberInpcrdFile(rst7_file)

    # Create the system
    system = prmtop.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.2*nanometers,
        constraints=None,
        implicitSolvent=app.OBC2
    )

    # Integrator with a 0.7 fs time step
    integrator = LangevinIntegrator(
        300*kelvin,      # Temperature
        1/picosecond,    # Friction coefficient
        0.7*femtoseconds # Time step
    )

    # Set up the simulation
    simulation = app.Simulation(prmtop.topology, system, integrator)
    simulation.context.setPositions(inpcrd.positions)

    # Minimize the energy
    simulation.minimizeEnergy()

    # Set up reporters
    simulation.reporters.append(app.DCDReporter(dcd_file, int(10*femtoseconds / 0.7*femtoseconds)))
    simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

    # Run the simulation for 100 ps
    n_steps = int(100*picoseconds / (0.7*femtoseconds))
    simulation.step(n_steps)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run OpenMM simulation with implicit solvent.')
    parser.add_argument('prmtop', type=str, help='Amber parameter/topology file (prmtop)')
    parser.add_argument('rst7', type=str, help='Amber restart file (rst7)')
    parser.add_argument('dcd', type=str, help='Output DCD trajectory file')
    args = parser.parse_args()

    run_simulation(args.prmtop, args.rst7, args.dcd)
