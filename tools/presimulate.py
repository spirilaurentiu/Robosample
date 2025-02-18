import argparse
from simtk.openmm import app, Platform, unit
from simtk.openmm import LangevinIntegrator, AndersenThermostat
from sys import stdout

def run_simulation(prmtop_file, rst7_file, dcd_file, simulation_time_hours):
    # Load the topology and positions
    prmtop = app.AmberPrmtopFile(prmtop_file)
    inpcrd = app.AmberInpcrdFile(rst7_file)

    # Create the system with HBonds constraints
    system = prmtop.createSystem(
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometers,
        constraints=app.HBonds,
        implicitSolvent=app.OBC2
    )

    # Add the Andersen thermostat to the system
    temperature = 300 * unit.kelvin
    collision_frequency = 1 / unit.picoseconds
    andersen_thermostat = AndersenThermostat(temperature, collision_frequency)
    system.addForce(andersen_thermostat)

    timestep = 2 * unit.femtoseconds
    equil_time = 5 * unit.nanoseconds
    write_time = 10 * unit.picoseconds

    integrator = LangevinIntegrator(
        temperature,        # Temperature
        1 / unit.picoseconds,  # Friction coefficient
        timestep            # Time step
    )

    # Set up the simulation
    platform = Platform.getPlatformByName('CUDA')
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    # Minimize the energy
    simulation.minimizeEnergy()

    # Equilibration run (5 ns)
    equilibration_steps = int(equil_time / timestep)
    simulation.step(equilibration_steps)

    # Remove previous reporters and add new ones for the main simulation
    write_interval = int(write_time / timestep)
    simulation.reporters.append(app.DCDReporter(dcd_file, write_interval))

    # Main simulation
    simulation.runForClockTime(simulation_time_hours * unit.hour)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run OpenMM simulation with implicit solvent.')
    parser.add_argument('prmtop', type=str, help='Amber parameter/topology file (prmtop)')
    parser.add_argument('rst7', type=str, help='Amber restart file (rst7)')
    parser.add_argument('dcd', type=str, help='Output DCD trajectory file')
    parser.add_argument('simulation_time_hours', type=float, help='Simulation time in hours')
    args = parser.parse_args()

    run_simulation(args.prmtop, args.rst7, args.dcd, args.simulation_time_hours)



# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1A5E.prmtop ../datasets/diverse-cath-nmr/1A5E_min.rst7 1A5E.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2JR6.prmtop ../datasets/diverse-cath-nmr/2JR6_min.rst7 2JR6.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2L29.prmtop ../datasets/diverse-cath-nmr/2L29_min.rst7 2L29.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1APQ.prmtop ../datasets/diverse-cath-nmr/1APQ_min.rst7 1APQ.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KB1.prmtop ../datasets/diverse-cath-nmr/2KB1_min.rst7 2KB1.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2L3B.prmtop ../datasets/diverse-cath-nmr/2L3B_min.rst7 2L3B.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1BLA.prmtop ../datasets/diverse-cath-nmr/1BLA_min.rst7 1BLA.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KCU.prmtop ../datasets/diverse-cath-nmr/2KCU_min.rst7 2KCU.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LGJ.prmtop ../datasets/diverse-cath-nmr/2LGJ_min.rst7 2LGJ.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/1L1I.prmtop ../datasets/diverse-cath-nmr/1L1I_min.rst7 1L1I.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KNQ.prmtop ../datasets/diverse-cath-nmr/2KNQ_min.rst7 2KNQ.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LT2.prmtop ../datasets/diverse-cath-nmr/2LT2_min.rst7 2LT2.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2EZA.prmtop ../datasets/diverse-cath-nmr/2EZA_min.rst7 2EZA.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KPU.prmtop ../datasets/diverse-cath-nmr/2KPU_min.rst7 2KPU.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2LTD.prmtop ../datasets/diverse-cath-nmr/2LTD_min.rst7 2LTD.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2HHI.prmtop ../datasets/diverse-cath-nmr/2HHI_min.rst7 2HHI.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KQ8.prmtop ../datasets/diverse-cath-nmr/2KQ8_min.rst7 2KQ8.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2N9B.prmtop ../datasets/diverse-cath-nmr/2N9B_min.rst7 2N9B.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2JOZ.prmtop ../datasets/diverse-cath-nmr/2JOZ_min.rst7 2JOZ.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/2KTA.prmtop ../datasets/diverse-cath-nmr/2KTA_min.rst7 2KTA.dcd 16 &
# nohup python3 presimulate.py ../datasets/diverse-cath-nmr/6V88.prmtop ../datasets/diverse-cath-nmr/6V88_min.rst7 6V88.dcd 16 &
