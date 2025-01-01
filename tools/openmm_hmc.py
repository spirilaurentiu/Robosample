import sys
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import HMCIntegrator
import argparse

def simulate(prmtop_file, rst7_file, dcd_file):
    # Load the prmtop and rst7 files
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(rst7_file)

    # Create a simulation system with 1.2 nm cutoff
    system = prmtop.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.2*nanometers,
        constraints = None,
        rigidWater = True,
        removeCMMotion = False,
        implicitSolvent=OBC2
    )

    # Set up a thermostat with stronger coupling
    system.addForce(AndersenThermostat(300 * kelvin, 5 / picosecond))

    # Set up GBSA with 1.2 nm cutoff
    gbsa = GBSAOBCForce()
    gbsa.setNonbondedMethod(GBSAOBCForce.CutoffNonPeriodic)
    gbsa.setCutoffDistance(1.2 * nanometer)
    gbsa.setSoluteDielectric(1.0)
    gbsa.setSolventDielectric(78.5)

    # Add particles to GBSA force
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded = force
            break
        if isinstance(force, GBSAOBCForce):
            force.setSoluteDielectric(1.0)
            force.setSolventDielectric(78.5)
    
    for i in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        radius = sigma * 0.5 * 2**(1/6.0)
        gbsa.addParticle(charge, radius, 1.0)

    system.addForce(gbsa)

    # Set up an integrator with smaller timestep
    temperature = 300 * kelvin
    step_size = 0.0007 * picoseconds
    steps_per_hmc_move = 1429

    # Define an HMC integrator
    integrator = HMCIntegrator(
        temperature=temperature,
        nsteps=steps_per_hmc_move,
        timestep=step_size
    )

    # Create a simulation object with CPU platform
    platform = Platform.getPlatformByName('CPU') # OpenCL
    simulation = Simulation(prmtop.topology, system, integrator, platform)

    # Set the positions from the rst7 file
    simulation.context.setPositions(inpcrd.positions)

    # If velocities are available, set them as well
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    # More thorough energy minimization
    print("Minimizing energy...")
    simulation.minimizeEnergy()  # Increased iterations
    print("Energy minimization complete.")

    # Initialize velocities at 10K
    simulation.context.setVelocitiesToTemperature(10*kelvin)

    # Add reporters for more frequent monitoring during heating
    simulation.reporters.append(StateDataReporter(sys.stdout, 50,
        step=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True  # Added speed monitoring
    ))

    # Gradual heating to 300K
    print("Heating system gradually...")
    for temp in [50, 100, 200, 300]:
        # Update thermostat temperature
        for force in system.getForces():
            if isinstance(force, AndersenThermostat):
                force.setDefaultTemperature(temp*kelvin)
                break
        
        # Run for 1000 steps at each temperature
        simulation.step(100)
        print(f"Completed heating to {temp}K")

    # Production run
    print("Running production simulation...")
    simulation.reporters.append(DCDReporter(dcd_file, 50))
    simulation.step(100000)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Simulate molecular dynamics.')
    parser.add_argument('--prmtop', required=True, help='Path to the prmtop file')
    parser.add_argument('--rst7', required=True, help='Path to the rst7 file')
    parser.add_argument('--dcd', required=True, help='Path to the dcd file')

    args = parser.parse_args()

    simulate(args.prmtop, args.rst7, args.dcd)


# python3 openmm_hmc.py --prmtop ../datasets/diverse-cath-nmr/1A5E.prmtop --rst7 ../datasets/diverse-cath-nmr/1A5E_min.rst7 --dcd 1A5E.dcd
# python3 openmm_hmc.py --prmtop ../datasets/alanine-dipeptide/alanine-dipeptide.prmtop --rst7 ../datasets/alanine-dipeptide/alanine-dipeptide.rst7 --dcd alanine-dipeptide.dcd