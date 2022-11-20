import numpy as np
import mdtraj as md
import flexor

prmtop = "../../tests_inputs/ligand.prmtop"
inpcrd = "../../tests_inputs/ligand.min.rst7"

## Load system
mdtrajObj = md.load(inpcrd, top=prmtop)

## Instantiate Flexor object
flexorObj = flexor.Flexor(mdtrajObj)

## Load DB file into Flexor object
flexorObj.loadFlexDB("./databases/aminoacids.flex.txt")

## Get flexibility file.
flexorObj.addWorld(range="resid 0 to 10", distanceCutoff=0, subset=["all"],
                 type="stretch" ,jointType="Pin")