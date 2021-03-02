# Robosample tools
# Works with Python 3
from __future__ import print_function, division
import mdtraj as md
import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry, distance, dihedral
import warnings


import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--top', type=str, 
help='topology file of system (prmtop)')
parser.add_argument('--traj', type=str, 
help='file containing coordinates (dcd, rst7)')
parser.add_argument('--probesize', type=float, 
help='Accesibility probe size in nm.')
parser.add_argument('--flex', type=str, 
help='file containing flexibilities, must correspond to prmtop/inpcrd selected')
parser.add_argument('--TD', action='store_true', default=False,
  help='')
args = parser.parse_args()

# Get
#lines = np.loadtxt(args.top)
infile = open(args.top, 'r')
lines = infile.readlines()
infile.close()

flag = False
dihlix = 0
dihIxs = []
for lix in range(len(lines)):
	line = lines[lix].rstrip().split()
	# Start recording when dihedrals are met
	if len(line) >= 2:
		if (line[0] == '%FLAG'):
			flag = False
		if line[1] == 'DIHEDRALS_INC_HYDROGEN':
			flag = True
	if flag == True:
		dihlix = dihlix + 1
		if dihlix > 2:
			# A = N/3 + 1
			onedihIxs = np.array(line[0:4], dtype=float)
			if onedihIxs[2] >= 0:
				onedihIxs = np.abs(onedihIxs)
				onedihIxs = (onedihIxs / 3.0).astype(int)
				dihIxs.append(onedihIxs)
dihIxs = np.array(dihIxs)

# Get molecule
molecule = md.load(args.traj, top=args.top)

# Get topology
topology = molecule.topology
table, bonds = topology.to_dataframe()

# Compute accesibility: shape (n_frames,n_atoms)
sasa = md.shrake_rupley(molecule, probe_radius = args.probesize, mode = 'residue') # nm
#for i in range(sasa.shape[1]):
#	print(i+1, sasa[0][i])
#exit(0)

#print("#i1 i2 Joint # name1 elem1 resid1 resname1 name2 elem2 resid2 resname2 sasa1 sasa2")
if(args.TD):
	for i in range(dihIxs.shape[0]):
		print("%5d %5d Pin # %4s %s %6d %3s   %4s %s %6d %3s %8.5f %8.5f" % (dihIxs[i][1], dihIxs[i][2] \
			,table.values[ dihIxs[i][1] ][1], table.values[ dihIxs[i][1] ][2] \
			,table.values[ dihIxs[i][1] ][3], table.values[ dihIxs[i][1] ][4] \
			,table.values[ dihIxs[i][2] ][1], table.values[ dihIxs[i][2] ][2] \
			,table.values[ dihIxs[i][2] ][3], table.values[ dihIxs[i][2] ][4] \
			,sasa[0][ dihIxs[i][1] ], sasa[0][ dihIxs[i][2] ]
		))
else:
	for i in range(topology.n_bonds):
		print("%5d %5d Cartesian # %4s %s %6d %3s   %4s %s %6d %3s %8.5f %8.5f" % (bonds[i][0], bonds[i][1] \
			,table.values[ int(bonds[i][0]) ][1], table.values[ int(bonds[i][0]) ][2] \
			,table.values[ int(bonds[i][0]) ][3], table.values[ int(bonds[i][0]) ][4] \
			,table.values[ int(bonds[i][1]) ][1], table.values[ int(bonds[i][1]) ][2] \
			,table.values[ int(bonds[i][1]) ][3], table.values[ int(bonds[i][1]) ][4] \
			,sasa[0][ table.values[ int(bonds[i][0]) ][3] ], sasa[0][ table.values[ int(bonds[i][1]) ][3] ]
		))


# Get vicinity
#query_indices = np.array(range(171))
#radius = 0.5 # nm
#neighbors = np.array(md.compute_neighbors(molecule, radius, query_indices))
#
#for i in range(neighbors.shape[0]):
#	for j in range(neighbors.shape[1]):
#		pass
#		#print(neighbors[i][j], end=' ')
#	#print(end='\n')

