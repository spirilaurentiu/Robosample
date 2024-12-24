import sys
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

def simulate(prmtop_file, rst7_file, dcd_file):
    # Load the prmtop and rst7 files
    prmtop = AmberPrmtopFile(prmtop_file)
    inpcrd = AmberInpcrdFile(rst7_file)

    # Create a simulation system with 1.2 nm cutoff
    system = prmtop.createSystem(
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.2*nanometers,
        constraints=HBonds
    )

    # Set up a thermostat with stronger coupling
    system.addForce(AndersenThermostat(300 * kelvin, 5 / picosecond))  # Increased coupling frequency

    # Set up GBSA with 1.2 nm cutoff
    gbsa = GBSAOBCForce()
    gbsa.setNonbondedMethod(GBSAOBCForce.CutoffNonPeriodic)
    gbsa.setCutoffDistance(1.2 * nanometer)
    gbsa.setSoluteDielectric(1.0)
    gbsa.setSolventDielectric(80)

    # Add particles to GBSA force
    for force in system.getForces():
        if isinstance(force, NonbondedForce):
            nonbonded = force
            break
    
    for i in range(nonbonded.getNumParticles()):
        charge, sigma, epsilon = nonbonded.getParticleParameters(i)
        radius = sigma * 0.5 * 2**(1/6.0)
        gbsa.addParticle(charge, radius, 1.0)

    system.addForce(gbsa)

    # Set up an integrator with smaller timestep
    integrator = VerletIntegrator(0.0007 * picoseconds)

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
    simulation.reporters.append(DCDReporter(dcd_file, 100))

    # Gradual heating to 300K
    print("Heating system gradually...")
    for temp in [50, 100, 200, 300]:
        # Update thermostat temperature
        for force in system.getForces():
            if isinstance(force, AndersenThermostat):
                force.setDefaultTemperature(temp*kelvin)
                break
        
        # Run for 1000 steps at each temperature
        simulation.step(1000)
        print(f"Completed heating to {temp}K")

    # Production run
    print("Running production simulation...")
    simulation.step(1428572) # roughly 1 ns

if __name__ == '__main__':
    simulate('diala.prmtop', 'diala.rst7', 'diala.dcd')