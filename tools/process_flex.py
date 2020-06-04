from ls_parsetxt import ParseTxt
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('--inFN', default=None, 
  help='Input file.')
parser.add_argument('--subset', default="phi", 
  help='Subset of bonds to extract.')
parser.add_argument('--accRange', type=float, default=[0.0, 1.0], nargs='+',
  help='Accesibility range.')
args = parser.parse_args()

# Read data
pars = ParseTxt()
pars.Read(args.inFN)
pdata = pars.parsed_data

aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

# Get phi
if (args.subset).lower() == "phi":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if(line[7] != "PRO"):
			if(((line[4] == "N") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "N"))):
				for word in line: print word,
				print

# Get psi
if (args.subset).lower() == "psi":
	for pi in range(len(pdata)):
		line = pdata[pi]
		if(line[7] != "PRO"):
			if(((line[4] == "C") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "C"))):
				for word in line: print word,
				print

# Get sidechains
if (args.subset).lower() == "sidechain":
	for pi in range(len(pdata)):
		line = pdata[pi]
		atom1Acc = float(line[-2])
		atom2Acc = float(line[-1])
		#print atom1Acc, args.accRange[0], args.accRange[1]
		if( ((atom1Acc >= args.accRange[0]) and (atom1Acc <= args.accRange[1])) or 
		    ((atom2Acc >= args.accRange[0]) and (atom2Acc <= args.accRange[1])) ):
			if(line[7] == "ALA"):
				if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
					for word in line: print word,
					print
			if(line[7] == "VAL"):
				if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
					for word in line: print word,
					print
				if(((line[4] == "CB") and (line[8] == "CG1")) or ((line[4] == "CG1") and (line[8] == "CB")) or 
				   ((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))):
					for word in line: print word,
					print
					



