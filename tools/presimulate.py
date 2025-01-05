import argparse
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
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
    platform = Platform.getPlatformByName('CUDA')
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    # Minimize the energy
    simulation.minimizeEnergy()

    # Set up reporters
    simulation.reporters.append(app.DCDReporter(dcd_file, int((10*femtoseconds) / (0.7*femtoseconds))))
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


# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1A5E.prmtop ../datasets/diverse-cath-nmr/1A5E_min.rst7 1A5E.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2JR6.prmtop ../datasets/diverse-cath-nmr/2JR6_min.rst7 2JR6.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2L29.prmtop ../datasets/diverse-cath-nmr/2L29_min.rst7 2L29.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1APQ.prmtop ../datasets/diverse-cath-nmr/1APQ_min.rst7 1APQ.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KB1.prmtop ../datasets/diverse-cath-nmr/2KB1_min.rst7 2KB1.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2L3B.prmtop ../datasets/diverse-cath-nmr/2L3B_min.rst7 2L3B.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1BLA.prmtop ../datasets/diverse-cath-nmr/1BLA_min.rst7 1BLA.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KCU.prmtop ../datasets/diverse-cath-nmr/2KCU_min.rst7 2KCU.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LGJ.prmtop ../datasets/diverse-cath-nmr/2LGJ_min.rst7 2LGJ.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1L1I.prmtop ../datasets/diverse-cath-nmr/1L1I_min.rst7 1L1I.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KNQ.prmtop ../datasets/diverse-cath-nmr/2KNQ_min.rst7 2KNQ.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LT2.prmtop ../datasets/diverse-cath-nmr/2LT2_min.rst7 2LT2.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2EZA.prmtop ../datasets/diverse-cath-nmr/2EZA_min.rst7 2EZA.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KPU.prmtop ../datasets/diverse-cath-nmr/2KPU_min.rst7 2KPU.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LTD.prmtop ../datasets/diverse-cath-nmr/2LTD_min.rst7 2LTD.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2HHI.prmtop ../datasets/diverse-cath-nmr/2HHI_min.rst7 2HHI.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KQ8.prmtop ../datasets/diverse-cath-nmr/2KQ8_min.rst7 2KQ8.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2N9B.prmtop ../datasets/diverse-cath-nmr/2N9B_min.rst7 2N9B.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2JOZ.prmtop ../datasets/diverse-cath-nmr/2JOZ_min.rst7 2JOZ.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KTA.prmtop ../datasets/diverse-cath-nmr/2KTA_min.rst7 2KTA.dcd &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/6V88.prmtop ../datasets/diverse-cath-nmr/6V88_min.rst7 6V88.dcd &
