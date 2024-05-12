## This is an example script for using the flexor2 module.

## Imports
import flexor2
import mdtraj as md

## Select the files for the system
prmtop = "../tests_inputs/2but/ligand.prmtop"
inpcrd = "../tests_inputs/2but/ligand.rst7"

## Instantiate the MDTraj Trajectory object
mdtrajObj = md.load(inpcrd, top=prmtop)

## Instantiate the Flexor object
flexorObj = flexor2.Flexor(mdtrajObj)

## Create the worlds

## The signature of the addWorld function is as follows:
#   range           =   DSL selection, identical to the one used by MDTraj
#   subset          =   list of strings, which can contain ("all", "side", "rama", "phi", "psi", "omega").
#                       This further filters the selections
#   jointType       =   "Cartesian", "Pin" or "Ball". Type of joint to use
# 
#   distanceCutoff  =   Include other atoms in the system that are within this distance
#                       from any atom in the selection defined by range. Defaults to 0.0
#   sasa_value      =   If supplied, will filter atoms from range so that only atoms that have a higher
#                       SASA than this are kept. Defaults to -1 (off)
#   rolling         =   If True, this world's flexibility will be split across multiple files (one joint per file)
#                       Defaults to False
#   FNOut           =   Root name for the resulting flex file(s)


## Simple examples of the function calls:

## Fully flexible world:
a = flexorObj.addWorld(range="all", subset=["all"],
                           jointType="Cartesian",
                           FNOut="./2but.all.cartesian")

## TD World:
a = flexorObj.addWorld(range="all", subset=["all"],
                           jointType="Pin",
                           FNOut="./2but.all.pin")

## TD backbone world: 
a = flexorObj.addWorld(range="all", subset=["rama"],
                           jointType="Pin",
                           FNOut="./2but.rama.pin")
