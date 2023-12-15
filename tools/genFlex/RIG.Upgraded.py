import numpy as np
import mdtraj as md
import flexor

roots = ["./chainA", "./chainB", "./chainC", "./chainD",
           "./chainE", "./chainF", "./chainG", "./chainH",
           "./DNA_chain1", "./DNA_chain2"]

roll_seles =    ["resid 1 to 45", "resid 1 to 21", "resid 1 to 16 or resid 117 to 130", "resid 1 to 38",
                "resid 1 to 45", "resid 1 to 26", "resid 1 to 17 or resid 117 to 130", "resid 1 to 38", 
                "none", "none"]
## We know apriori how many molecules we'll have
allWorlds_flex = [[] for x in range(len(roots))] 

for moleculeIx in range(len(roots)):
## Load system
    moleculeWorlds = []
    mdtrajObj = md.load("{}.rst7".format(roots[moleculeIx]), 
                        top="{}.prmtop".format(roots[moleculeIx]))

    ## Instantiate Flexor object
    flexorObj = flexor.Flexor(mdtrajObj)

    ## Generate Cartesian Flex Files (all DoFs)
    a = flexorObj.addWorld(range="all", subset=["all"],
                           jointType="Cartesian", 
                           FNOut="./flexfiles/Nucleosome.{}.Cartesian".format(moleculeIx))
    moleculeWorlds.append(a)
    
    ## Generate Rama Pin Flex File 
    a = flexorObj.addWorld(range="all", subset=["rama"], jointType="Pin",
                           FNOut="./flexfiles/Nucleosome.{}.Rama.Pin".format(moleculeIx), rolling=False)
    moleculeWorlds.append(a)

    ## Generate Sidechain Pin Flex File 
    a = flexorObj.addWorld(range="all", subset=["side"], jointType="Pin",
                           FNOut="./flexfiles/Nucleosome.{}.Side.Pin".format(moleculeIx), rolling=False)
    moleculeWorlds.append(a)

    ## Generate Histone Tails Flex Files (Pin Rolls)
    a = flexorObj.addWorld(range=roll_seles[moleculeIx], subset=["rama"], jointType="Pin",
                           FNOut="./flexfiles/Nucleosome.{}.Tail.Pin".format(moleculeIx), rolling=True)
    moleculeWorlds.append(a)

    for dim1 in moleculeWorlds:
        for dim2 in dim1:
            allWorlds_flex[moleculeIx].append(dim2)

maxWorlds = 0
for _ in allWorlds_flex:
    if (len(_) > maxWorlds):
        maxWorlds = len(_)
#    print ("\n\n")
#    print (_)

## We add rigid flexibility files to have the same number of flex files
## for all molecules
for _ in allWorlds_flex:
    if (len(_) < maxWorlds):
        [_.append("rigid.flex") for x in range(maxWorlds - len(_))]

#for _ in allWorlds_flex:
#    print ("\n")
#    print (_)

## Now that we have the molecules and the associated flex files,
## we can generate the final input file
fileOut = open ("test.inp", "w")
fileOut.write("PRMTOP ")
for worldIx in range(maxWorlds):
    for _ in (roots):
        fileOut.write("{}.prmtop ".format(_))

fileOut.write("\nINPCRD ")
for worldIx in range(maxWorlds):
    for _ in (roots):
        fileOut.write("{} ".format(_))

fileOut.write("\nRBFILE ")
for _ in range(len(roots)*(maxWorlds)):
    fileOut.write("ligand.rb ")

fileOut.write("\nROOTS ")
for _ in range(len(roots)*(maxWorlds)):
    fileOut.write("0 ")

fileOut.write("\nROOT_MOBILITY ")
for _ in range(len(roots)*(maxWorlds)):
    fileOut.write("Weld ")

fileOut.write("\nOUTPUT_DIR temp")

fileOut.write("\nROUNTS_TILL_REBLOCK ")
for _ in range(maxWorlds):
    fileOut.write("10 ")

fileOut.write("\nWORLDS ")
for _ in range(maxWorlds):
    fileOut.write("R{} ".format(_))

fileOut.write("\nRANDOM_WORLD_ORDER ")
for _ in range(maxWorlds):
    fileOut.write("FALSE ")

fileOut.write("\nREPRODUCIBLE ")
for _ in range(maxWorlds):
    fileOut.write("FALSE ")

fileOut.write("\nDISTORT_OPTION ")
for _ in range(maxWorlds):
    fileOut.write("0 ")

fileOut.write("\nSAMPLES_PER_ROUND ")
for _ in range(maxWorlds):
    fileOut.write("1 ")

fileOut.write("\nFIXMAN_POTENTIAL ")
for _ in range(maxWorlds):
    fileOut.write("TRUE ")

fileOut.write("\nFIXMAN_TORQUE ")
for _ in range(maxWorlds):
    fileOut.write("TRUE ")

fileOut.write("\nGBSA ")
for _ in range(maxWorlds):
    fileOut.write("1.0 ")

fileOut.write("\nPRINT_FREQ ")
for _ in range(maxWorlds):
    fileOut.write("1 ")

fileOut.write("\nWRITEPDBS 1")
for _ in range(maxWorlds-1):
    fileOut.write("0 ")

fileOut.write("\nVISUAL ")
for _ in range(maxWorlds):
    fileOut.write("FALSE ")

fileOut.write("\nGEOMETRY ")
for _ in range(maxWorlds):
    fileOut.write("FALSE ")

fileOut.write("\nTHREADS ")
for _ in range(maxWorlds):
    fileOut.write("0 ")

fileOut.write("\nOPENMM ")
for _ in range(maxWorlds):
    fileOut.write("TRUE ")

fileOut.write("\nOPENMM_CalcOnlyNonbonded ")
for _ in range(maxWorlds):
    fileOut.write("FALSE ")

fileOut.write("\nNONBONDED_METHOD ")
for _ in range(maxWorlds):
    fileOut.write("0 ")

fileOut.write("\nNONBONDED_CUTOFF ")
for _ in range(maxWorlds):
    fileOut.write("1.2 ")

fileOut.write("\nMOLECULES ")
for _ in range(len(roots)):
    fileOut.write("nucleosome ")

## The big one
fileOut.write("\nFLEXFILE ")
for worldIx in range(maxWorlds):
    for moleculeIx in range(len(roots)):
        fileOut.write("{} ".format(allWorlds_flex[moleculeIx][worldIx]))

fileOut.write("\nRUN_TYPE ")
for _ in range(maxWorlds):
    fileOut.write("NORMAL ")

fileOut.write("\nFLOW_OPTION ")
for _ in range(maxWorlds):
    fileOut.write("0 ")

fileOut.write("\nWORK_OPTION ")
for _ in range(maxWorlds):
    fileOut.write("0 ")

fileOut.write("\nROUNDS 10000 ")

fileOut.write("\nSEED ")
for _ in range(maxWorlds):
    fileOut.write("451 ")

fileOut.write("\nTIMESTEPS ")
for _ in range(maxWorlds):
    fileOut.write("0.01 ")

fileOut.write("\nROUNDS_TILL_REBLOCK ")
for _ in range(maxWorlds):
    fileOut.write("1 ")

fileOut.write("\nMDSTEPS ")
for _ in range(maxWorlds):
    fileOut.write("1 ")

fileOut.write("\nBOOST_MDSTEPS ")
for _ in range(maxWorlds):
    fileOut.write("1 ")

fileOut.write("\nTEMPERATURE_INI ")
for _ in range(maxWorlds):
    fileOut.write("300 ")

fileOut.write("\nTEMPERATURE_FIN ")
for _ in range(maxWorlds):
    fileOut.write("300 ")

fileOut.write("\nBOOST_TEMPERATURE ")
for _ in range(maxWorlds):
    fileOut.write("300 ")

fileOut.write("\nTHERMOSTAT ")
for _ in range(maxWorlds):
    fileOut.write("Andersen ")

fileOut.write("\nFFSCALE ")
for _ in range(maxWorlds):
    fileOut.write("AMBER ")

fileOut.write("\nSAMPLERS ")
for _ in range(maxWorlds):
    fileOut.write("MC ")

fileOut.write("\nINTEGRATORS ")
for _ in range(maxWorlds):
    fileOut.write("VV ")


