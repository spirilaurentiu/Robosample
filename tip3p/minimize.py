from sys import stdout

from openmm.app import *
from openmm.unit import *

from openmm import *

# -------------------------
# Load AMBER system
# -------------------------
prmtop = AmberPrmtopFile("2ala.prmtop")
inpcrd = AmberInpcrdFile("2ala.rst7")

# -------------------------
# Build system (NVT first)
# -------------------------
system = prmtop.createSystem(
    nonbondedMethod=PME, nonbondedCutoff=10 * angstrom, constraints=HBonds
)
platform = Platform.getPlatformByName("CUDA")


def make_integrator():
    return LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 2 * femtoseconds)


# -------------------------
# MINIMIZATION + NVT
# -------------------------
NVT_STEPS = 250000  # 0.5 ns @ 2 fs

integrator = make_integrator()
simulation = Simulation(prmtop.topology, system, integrator, platform=platform)
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# live progress -> stdout
simulation.reporters.append(
    StateDataReporter(
        stdout,
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=NVT_STEPS,
        separator="\t",
    )
)

print("Minimizing...")
e0 = simulation.context.getState(
    getEnergy=True, enforcePeriodicBox=True
).getPotentialEnergy()
print("  PE before:", e0)
simulation.minimizeEnergy(maxIterations=5000)
e1 = simulation.context.getState(
    getEnergy=True, enforcePeriodicBox=True
).getPotentialEnergy()
print("  PE after: ", e1)

print("NVT equilibration...")
simulation.context.setVelocitiesToTemperature(300 * kelvin)
simulation.step(NVT_STEPS)

# Save NVT state (box vectors come for free, no kwarg needed)
state = simulation.context.getState(
    getPositions=True, getVelocities=True, enforcePeriodicBox=True
)
positions = state.getPositions()
velocities = state.getVelocities()
box = state.getPeriodicBoxVectors()

# -------------------------
# NPT SETUP (NEW CONTEXT)
# -------------------------
NPT_STEPS = 1000000  # 2 ns @ 2 fs

print("Switching to NPT...")
system.addForce(MonteCarloBarostat(1 * bar, 300 * kelvin, 25))

integrator = make_integrator()  # a context needs its own integrator
simulation = Simulation(prmtop.topology, system, integrator, platform=platform)
simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)
simulation.context.setPeriodicBoxVectors(*box)

simulation.reporters.append(
    StateDataReporter(
        stdout,
        5000,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        progress=True,
        remainingTime=True,
        speed=True,
        totalSteps=NPT_STEPS,
        separator="\t",
    )
)

print("NPT equilibration...")
simulation.step(NPT_STEPS)

# -------------------------
# FINAL STATE
# -------------------------
state = simulation.context.getState(
    getPositions=True, getVelocities=True, enforcePeriodicBox=True
)
positions = state.getPositions()
velocities = state.getVelocities()
box_vectors = state.getPeriodicBoxVectors()

a = box_vectors[0][0].value_in_unit(nanometer)
b = box_vectors[1][1].value_in_unit(nanometer)
c = box_vectors[2][2].value_in_unit(nanometer)
print("Final box size (nm):", a, b, c)

# -------------------------
# SAVE STATE
# -------------------------
import parmed

parm = parmed.load_file("2ala.prmtop", "2ala.rst7")
parm.positions = positions
parm.box_vectors = box_vectors
parm.save("equilibrated.rst7", overwrite=True)

print("Done. Wrote equilibrated_state.xml and equilibrated.rst7")
