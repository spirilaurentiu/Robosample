# Robosample tools
from ls_parsetxt import ParseTxt
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFN', default=None, 
  help='Input file.')
parser.add_argument('--subset', default="phi", 
  help='Subset of bonds to extract.')
parser.add_argument('--residRange', type=int, default=[0, 10000000], nargs='+',
  help='Residue index range.')
parser.add_argument('--accRange', type=float, default=[0.0, 1.0], nargs='+',
  help='Accesibility range.')
parser.add_argument('--joint', type=str, default='Pin',
  help='Joint type.')
args = parser.parse_args()

# Read data
pars = ParseTxt()
pars.Read(args.inFN)
pdata = pars.parsed_data

def printFlexLine(line, joint):
	print ("%d %d %s %s %s %s %d %s %s %s %d %s %.5f %.5f" % \
(int(line[0]), int(line[1]), joint, line[3], line[4], line[5], int(line[6]), \
line[7], line[8], line[9], int(line[10]), line[11], float(line[12]), \
float(line[13]))
	)
#

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
sugars = ['UYB', '4YB', 'VMB', '2MA', '3LB', '0fA', '0SA']
# Get all
if (args.subset).lower() == "all":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				printFlexLine(line, args.joint)

# Get phi
if (args.subset).lower() == "phi":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				if(line[7] != "PRO"):
					if(((line[4] == "N") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "N"))):
						printFlexLine(line, args.joint)

# Get psi
if (args.subset).lower() == "psi":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				if(((line[4] == "C") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "C"))):
					printFlexLine(line, args.joint)

# Get Ramachandran flexibility
if (args.subset).lower() == "rama":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				if(line[7] != "PRO"):
					if(((line[4] == "N") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "N"))):
						printFlexLine(line, args.joint)
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				if(((line[4] == "C") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "C"))):
					printFlexLine(line, args.joint)


# Get sidechains
if (args.subset).lower() == "side":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			line = pdata[pi]
			atom1Acc = float(line[-2])
			atom2Acc = float(line[-1])
			if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
			    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
				if(line[7] == "ALA"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
				if(line[7] == "VAL"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG1")) or ((line[4] == "CG1") and (line[8] == "CB")) or 
					   ((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))):
						printFlexLine(line, args.joint)
				if(line[7] == "LEU"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "CD1")) or ((line[4] == "CD1") and (line[8] == "CG")) or 
					   ((line[4] == "CG") and (line[8] == "CD2")) or ((line[4] == "CD2") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
				if(line[7] == "ILE"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG1")) or ((line[4] == "CG1") and (line[8] == "CB")) or 
					   ((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CG1") and (line[8] == "CD1")) or ((line[4] == "CD1") and (line[8] == "CG1"))):
						printFlexLine(line, args.joint)
				if(line[7] == "MET"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)

					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "CG") and (line[8] == "SD")) or ((line[4] == "SD") and (line[8] == "CG"))):
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "CG") and (line[8] == "SD")) or ((line[4] == "SD") and (line[8] == "CG"))):
							printFlexLine(line, args.joint)

					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "SD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "SD"))):
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "SD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "SD"))):
							printFlexLine(line, args.joint)

				if(line[7] == "PHE"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
				if(line[7] == "TRP"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
				if(line[7] == "TYR"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)

					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "CZ") and (line[8] == "OH")) or ((line[4] == "OH") and (line[8] == "CZ"))): 
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "CZ") and (line[8] == "OH")) or ((line[4] == "OH") and (line[8] == "CZ"))): 
							printFlexLine(line, args.joint)

				if(line[7] == "ASN"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
				if(line[7] == "SER"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)

					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "CB") and (line[8] == "OG")) or ((line[4] == "OG") and (line[8] == "CB"))): 
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "CB") and (line[8] == "OG")) or ((line[4] == "OG") and (line[8] == "CB"))): 
							printFlexLine(line, args.joint)
				if(line[7] == "THR"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "CB") and (line[8] == "OG1")) or ((line[4] == "OG1") and (line[8] == "CB"))):
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "CB") and (line[8] == "OG1")) or ((line[4] == "OG1") and (line[8] == "CB"))):
							printFlexLine(line, args.joint)
				if(line[7] == "CYS"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if args.joint in ["BallM", "BallF"] :
						if(((line[4] == "CB") and (line[8] == "SG")) or ((line[4] == "SG") and (line[8] == "CB"))): 
							printFlexLine(line, "Pin")
					else:
						if(((line[4] == "CB") and (line[8] == "SG")) or ((line[4] == "SG") and (line[8] == "CB"))): 
							printFlexLine(line, args.joint)
				if(line[7] == "GLN"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
				if(line[7] == "LYS"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "CD"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CE") and (line[8] == "NZ")) or ((line[4] == "NZ") and (line[8] == "CE"))):
						printFlexLine(line, args.joint)
				if(line[7] == "ARG"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CD") and (line[8] == "NE")) or ((line[4] == "NE") and (line[8] == "CD"))):
						printFlexLine(line, args.joint)
					#if(((line[4] == "NE") and (line[8] == "CZ")) or ((line[4] == "CZ") and (line[8] == "NE"))):
					#	printFlexLine(line, args.joint)
					#if(((line[4] == "CZ") and (line[8] == "NH1")) or ((line[4] == "NH1") and (line[8] == "CZ"))):
					#	printFlexLine(line, args.joint)
					#if(((line[4] == "CZ") and (line[8] == "NH2")) or ((line[4] == "NH2") and (line[8] == "CZ"))):
					#	printFlexLine(line, args.joint)
				if((line[7][0] == "H") and (line[7][1] == "I")):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
				if(line[7] == "ASP"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
				if(line[7] == "GLU"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
						
	
# Glycan part	
if (args.subset).lower() == "sugnln":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			if(line[7] == "NLN") and (line[11] == "NLN"):
					if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
						printFlexLine(line, args.joint)
					if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
						printFlexLine(line, args.joint)
					if(((line[4] == "CG") and (line[8] == "ND2")) or ((line[4] == "ND2") and (line[8] == "CG"))):
						printFlexLine(line, args.joint)
			if(line[7] == "NLN") and (line[11] in sugars):
				if((line[4] == "ND2") or (line[8] == "ND2")):
					printFlexLine(line, args.joint)
	
if (args.subset).lower() == "suginter":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			if(line[7] in sugars) and (line[11] in sugars) and (line[7] != line[11]):
				printFlexLine(line, args.joint)

if (args.subset).lower() == "sugout":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if((int(line[6]) >= args.residRange[0]) and (int(line[6]) <= args.residRange[1])):
			if(line[7] in sugars) and (line[7] == line[11]): # within the same sugar
				if(line[7] not in ['0SA']):
					if(((line[4] in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']) and (line[8] not in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'])) or \
					   ((line[8] in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']) and (line[4] not in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']))):
						if((line[5] != 'H') and (line[9] != 'H')):
							printFlexLine(line, args.joint)
				elif(line[7] in ['0SA']):
					if(((line[4] in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']) and (line[8] not in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6'])) or \
					   ((line[8] in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']) and (line[4] not in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']))):
						if((line[5] != 'H') and (line[9] != 'H')):
							printFlexLine(line, args.joint)






