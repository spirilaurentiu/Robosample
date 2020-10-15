# Works with Python 3
from __future__ import print_function, division

import sys, os, glob
import numpy as np
import scipy
import scipy.stats
import argparse
from autocorFuncs import *

import mdtraj as md

from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry, distance, dihedral
import warnings

class TrajectoryAnalyzer:
	def __init__(self, topFN, molName, FNSeeds, simDirs):
		'''Initialize'''
		# Input parameters variables
		self.molName = molName
		self.topFN = topFN
		self.simDirs = simDirs
		self.FNSeeds = FNSeeds

		self.trajectories = [None] * len(FNSeeds)

		self.rmsds = [None] * len(FNSeeds)
		self.RGs = [None] * len(FNSeeds)
		self.SASAs = [None] * len(FNSeeds)
		self.totSASAs = [None] * len(FNSeeds)
		self.helicities1 = [None] * len(FNSeeds)
	#

	def ReadPdbs(self, verbose = True):
		'''
		Reads trajectories from pdb files
		'''
		# Main loop: iterates through seeds = simulation reps
		for seedi in range(len(self.FNSeeds)):
			# Get file names
			FNList = []
			for di in range(len(self.simDirs)):
				FNList.append(glob.glob(self.simDirs[di] + "pdbs/sb." + self.molName + str(self.FNSeeds[seedi]) + '.0.00*pdb'))
			FNList = list(np.array(FNList).flat)
			if verbose == True:
				print(FNList)
			nfiles = len(FNList)

			self.trajectories[seedi] = md.load(FNList, top = self.topFN, stride = 1)
	#

	def Pdbs2dcd(self):
		'''
		Generates a dcd file from the set of pdb files in seed directory
		TODO
		'''
		# Main loop: iterates through seeds = simulation reps
		for seedi in range(len(self.FNSeeds)):
			# Get file names
			FNList = []
			for di in range(len(self.simDirs)):
				FNList.append(glob.glob(self.simDirs[di] + "pdbs/sb." + self.molName + str(self.FNSeeds[seedi]) + '.0.00*pdb'))
			FNList = list(np.array(FNList).flat)
			if verbose == True:
				print(FNList)
			nfiles = len(FNList)
	#
	

	def ReadDcds(self, verbose = True):
		'''
		Reads trajectories form dcd files
		'''
		# Main loop: iterates through seeds = simulation reps
		for seedi in range(len(self.FNSeeds)):
			# Get file names
			FNList = []
			for di in range(len(self.simDirs)):
				FNList.append(glob.glob(self.simDirs[di] + "traj." + self.molName + str(self.FNSeeds[seedi]) + '.dcd')[0])
			FNList = list(np.array(FNList).flat)
			if verbose == True:
				print(FNList)
			nfiles = len(FNList)
			self.trajectories[seedi] = md.load(FNList, top = self.topFN, stride = 1)
	#

	def RMSD(self):
		for seedi in range(len(self.FNSeeds)):
			self.rmsds[seedi] = md.rmsd(self.trajectories[seedi], self.trajectories[seedi], 0)
			print("RMSD:")
			print(self.rmsds[seedi])
	#

	def RG(self):
		for seedi in range(len(self.FNSeeds)):
			self.RGs[seedi] = md.compute_rg(self.trajectories[seedi], masses = None)
			print("Radius of gyration:")
			print((self.RGs[seedi]).shape)
	#

	def SASA(self):
		for seedi in range(len(self.FNSeeds)):
			self.SASAs[seedi] = md.shrake_rupley(self.trajectories[seedi], probe_radius = 0.14, mode='residue', get_mapping=False)
			nframes = (self.SASAs[seedi]).shape[0]
			nres = (self.SASAs[seedi]).shape[1]
			self.totSASAs[seedi] = np.zeros(( nframes ))
			for i in range(nframes):
				self.totSASAs[seedi][i] = np.sum(self.SASAs[seedi][i])
			print("SASA:")
			print((self.totSASAs[seedi]).shape)
			#print(self.SASAs[seedi])
	#

	def Helicity(self):
		for seedi in range(len(self.FNSeeds)):
			dssps = md.compute_dssp(self.trajectories[seedi], simplified=True)
			self.helicities1[seedi] = np.zeros(( len(dssps) ))

			#for i in range(len(dssps)):
			#	print(dssps[i])

			for i in range(len(dssps)):
				H = 0
				for j in range(len(dssps[i])):
					if dssps[i][j] == "H": H += 1
				self.helicities1[seedi][i] = H / len(dssps[i])

			#print(self.helicities1[seedi])
			
	#
#


## Get
##lines = np.loadtxt(args.top)
#infile = open(args.top, 'r')
#lines = infile.readlines()
#infile.close()
#
#flag = False
#dihlix = 0
#dihIxs = []
#for lix in range(len(lines)):
#	line = lines[lix].rstrip().split()
#	# Start recording when dihedrals are met
#	if len(line) >= 2:
#		if (line[0] == '%FLAG'):
#			flag = False
#		if line[1] == 'DIHEDRALS_INC_HYDROGEN':
#			flag = True
#	if flag == True:
#		dihlix = dihlix + 1
#		if dihlix > 2:
#			# A = N/3 + 1
#			onedihIxs = np.array(line[0:4], dtype=float)
#			if onedihIxs[2] >= 0:
#				onedihIxs = np.abs(onedihIxs)
#				onedihIxs = (onedihIxs / 3.0).astype(int)
#				dihIxs.append(onedihIxs)
#dihIxs = np.array(dihIxs)
#
## Get molecule
#molecule = md.load(args.traj, top=args.top)
#topology = molecule.topology
#table, bonds = topology.to_dataframe()
#
## Compute accesibility: shape (n_frames,n_atoms)
#sasa = md.shrake_rupley(molecule, probe_radius = args.probesize, mode = 'residue') # nm
##for i in range(sasa.shape[1]):
##	print(i+1, sasa[0][i])
##exit(0)
#
##print("#i1 i2 Joint # name1 elem1 resid1 resname1 name2 elem2 resid2 resname2 sasa1 sasa2")
#if(args.TD):
#	for i in range(dihIxs.shape[0]):
#		print("%5d %5d Pin # %4s %s %6d %3s   %4s %s %6d %3s %8.5f %8.5f" % (dihIxs[i][1], dihIxs[i][2] \
#			,table.values[ dihIxs[i][1] ][1], table.values[ dihIxs[i][1] ][2] \
#			,table.values[ dihIxs[i][1] ][3], table.values[ dihIxs[i][1] ][4] \
#			,table.values[ dihIxs[i][2] ][1], table.values[ dihIxs[i][2] ][2] \
#			,table.values[ dihIxs[i][2] ][3], table.values[ dihIxs[i][2] ][4] \
#			,sasa[0][ dihIxs[i][1] ], sasa[0][ dihIxs[i][2] ]
#		))
#else:
#	for i in range(topology.n_bonds):
#		print("%5d %5d Cartesian # %4s %s %6d %3s   %4s %s %6d %3s %8.5f %8.5f" % (bonds[i][0], bonds[i][1] \
#			,table.values[ int(bonds[i][0]) ][1], table.values[ int(bonds[i][0]) ][2] \
#			,table.values[ int(bonds[i][0]) ][3], table.values[ int(bonds[i][0]) ][4] \
#			,table.values[ int(bonds[i][1]) ][1], table.values[ int(bonds[i][1]) ][2] \
#			,table.values[ int(bonds[i][1]) ][3], table.values[ int(bonds[i][1]) ][4] \
#			,sasa[0][ table.values[ int(bonds[i][0]) ][3] ], sasa[0][ table.values[ int(bonds[i][1]) ][3] ]
#		))
