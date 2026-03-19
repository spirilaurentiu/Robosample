import parmed as pmd
from openmm import app
import openmm as mm
from openmm import unit

prmtop_file = "examples/pro.prmtop"
rst7_file = "examples/pro.rst7"
out_rst7 = "examples/pro-min.rst7"

structure = pmd.load_file(prmtop_file, rst7_file)
system = structure.createSystem(
    nonbondedMethod=app.NoCutoff,
    nonbondedCutoff=1.0 * unit.nanometer,
    constraints=app.HBonds
)

integrator = mm.LangevinIntegrator(
    300 * unit.kelvin,
    1.0 / unit.picosecond,
    0.002 * unit.picoseconds
)

platform = mm.Platform.getPlatformByName("CUDA")

simulation = app.Simulation(
    structure.topology,
    system,
    integrator,
    platform
)

simulation.context.setPositions(structure.positions)
simulation.minimizeEnergy(maxIterations=1000)
state = simulation.context.getState(getPositions=True)

structure.positions = state.getPositions()
structure.save(out_rst7, format='rst7')
