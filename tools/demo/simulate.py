from robosample import *

# Load Amber files
prmtop = AmberPrmtopFile("ala10/ligand.prmtop")
inpcrd = AmberInpcrdFile("ala10/ligand.rst7")

# Create a Robosample system by calling createSystem on prmtop
system = prmtop.createSystem(createDirs = False,
	nonbondedMethod = "NoCutoff",
 	nonbondedCutoff = 1.0*nanometer,
 	constraints = None,
 	rigidWater = True,
 	implicitSolvent = True,
 	soluteDielectric = 1.0,
 	solventDielectric = 78.5,
 	removeCMMotion = False
)

integrator = HMCIntegrator(300*kelvin,   # Temperature of head bath
                           0.006*picoseconds) # Time step

platform = Platform.getPlatformByName('GPU')

simulation = Simulation(prmtop.topology, system, integrator, platform)

simulation.reporters.append(PDBReporter('robots/', 10))

simulation.context.setPositions(inpcrd.positions)

# run simulation
simulation.step(50)

