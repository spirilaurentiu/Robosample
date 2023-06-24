from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout


from simtk.openmm.app.internal import amber_file_parser

filename = sys.argv[1]
prmtop = amber_file_parser.PrmtopLoader(filename)


atomsnames = prmtop.getAtomNames()
atomstypes = prmtop.getAtomTypes()
NonbondTerms = prmtop.getNonbondTerms()
Interactions14 = prmtop.get14Interactions()
ExcludedPairs = prmtop.getExcludedAtoms()


print(">>> TEST amber_file_parser OPENMM >>> \n\n")

print(">>> Atoms >>> \n")

print("ind", "Atom", "type", "rvdw(nm)", "eps(KJ/mol)", sep="\t")
for k in range(len(atomsnames)):
    print(k, atomsnames[k], atomstypes[k], NonbondTerms[k][0], NonbondTerms[k][1], sep="\t")



print("\n>>> Interactions 1-4 >>> \n")

print("ind", "Atom1", "Atom1","chargeProduct", "rMin(nm)", "eps(KJ/mol)", sep="\t")
for k in range(len(Interactions14)):
    print(k, Interactions14[k][0], Interactions14[k][1], Interactions14[k][2], Interactions14[k][3], Interactions14[k][4], sep="\t")



print("\n>>> Excluded Pairs List >>> \n")


count=1;
print("ind", "Atom1", "Atom1", sep="\t")
for k in range( len(ExcludedPairs) ):
    for j in range( len(ExcludedPairs[k]) ):
        print(count, k, ExcludedPairs[k][j], sep="\t")
        count = count + 1

print("\n\n\n")





# prmtop = AmberPrmtopFile('ligand.prmtop')
# inpcrd = AmberInpcrdFile('ligand.inpcrd')


# system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
# integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
# simulation = Simulation(prmtop.topology, system, integrator)
# simulation.context.setPositions(inpcrd.positions)
# if inpcrd.boxVectors is not None:
#     simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
# simulation.minimizeEnergy()
# simulation.reporters.append(PDBReporter('2but/output.pdb', 1000))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# simulation.step(10000)
#
