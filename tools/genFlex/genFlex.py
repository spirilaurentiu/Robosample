import numpy as np
import mdtraj as md
import flexor

prmtop = "example.0.prmtop"
inpcrd = "example.0.rst7"

## Load system
mdtrajObj = md.load(inpcrd, top=prmtop)

## Instantiate Flexor object
flexorObj = flexor.Flexor(mdtrajObj)

## Get flexibility file.
flexorObj.addWorld(range="all", distanceCutoff=0, subset=["all"],
                 jointType="Cartesian", sasa_value=-1.0,
                 FNOut="./example0.cart", rolling=False)

flexorObj.addWorld(range="all", distanceCutoff=0, subset=["rama"],
                 jointType="Pin", sasa_value=-1.0,
                 FNOut="./example0.pin.rama", rolling=False)
