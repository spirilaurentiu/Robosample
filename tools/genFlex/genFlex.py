import numpy as np
import mdtraj as md
import flexor

prmtop = "EXAMPLE.prmtop"
inpcrd = "EXAMPLE.rst7"

## Load system
mdtrajObj = md.load(inpcrd, top=prmtop)

## Instantiate Flexor object
flexorObj = flexor.Flexor(mdtrajObj)

## Load DB file into Flexor object
flexorObj.loadFlexDB("./aminoacids.flex.txt")

## Get flexibility file.
flexorObj.addWorld(range="all", distanceCutoff=0, subset=["all"],
                 jointType="Cartesian", sasa_value=-1.0,
                 FNOut="./protein.cart.new", rolling=False)

