# Robosample tools
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
		if isinstance(self.FNSeeds, list) != True:
			self.FNSeeds = [self.FNSeeds]
		print('TrajectoryAnalyzer init FNSeeds', self.FNSeeds)

		self.trajectories = [None] * len(FNSeeds)

		self.distances = [None] * len(FNSeeds)
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
				print('TrajectoryAnalyzer simDirs', self.simDirs[di])
				print('TrajectoryAnalyzer ReadDcds() FNSeeds', self.FNSeeds)
				print('TrajectoryAnalyzer simDirs get', self.simDirs[di] + "traj." + self.molName + "." + str(self.FNSeeds[seedi]) + '.dcd')
				fnlist = glob.glob(self.simDirs[di] + "traj." + self.molName + "." + str(self.FNSeeds[seedi]) + '.dcd')
				print('TrajectoryAnalyzer simDirs got', fnlist, '. Taking the first one.')
				FNList.append(fnlist[0])
			FNList = list(np.array(FNList).flat)
			if verbose == True:
				print(FNList)
			nfiles = len(FNList)
			self.trajectories[seedi] = md.load(FNList, top = self.topFN, stride = 1)
	#

	def Distance(self, indeces):
		for seedi in range(len(self.FNSeeds)):
			self.distances[seedi] = md.compute_distances(self.trajectories[seedi], indeces)
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


