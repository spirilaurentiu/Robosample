from openmm import Context, Platform, VerletIntegrator, unit
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, NoCutoff

prmtop = AmberPrmtopFile("2but/ligand.prmtop")
inpcrd = AmberInpcrdFile("2but/ligand.rst7")

system = prmtop.createSystem(
    nonbondedMethod=NoCutoff,
    constraints=None,
    rigidWater=False,
)

# Give each force a separate group for diagnosis.
for i, force in enumerate(system.getForces()):
    force.setForceGroup(i)

integrator = VerletIntegrator(0.001 * unit.picoseconds)
context = Context(system, integrator, Platform.getPlatformByName("Reference"))
context.setPositions(inpcrd.positions)

state = context.getState(getEnergy=True)
print("Total:", state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole))

for i, force in enumerate(system.getForces()):
    energy = context.getState(getEnergy=True, groups={i}).getPotentialEnergy()
    print(i, type(force).__name__, energy.value_in_unit(unit.kilojoules_per_mole))
