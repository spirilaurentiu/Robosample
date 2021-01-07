from robosample import *

# Load Amber files
prmtop = AmberPrmtopFile("ala10/ligand.prmtop")
inpcrd = AmberInpcrdFile("ala10/ligand.rst7")

# Create a Robosample system by calling createSystem on prmtop
system = prmtop.createSystem(createDirs = False)
integrator = HMCIntegrator(300*kelvin,   # Temperature of head bath
                           0.002*picoseconds) # Time step

simulation = Simulation(prmtop.topology, system, integrator)
simulation.context.setPositions(inpcrd.positions)

# run simulation
simulation.step(50)

