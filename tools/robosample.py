# Robosample tools
# Imports
import os, sys, errno, subprocess
import shutil
import glob
import itertools

import mdtraj as md
import numpy as np

from mdtraj.utils import ensure_type
from mdtraj.geometry import _geometry, distance, dihedral

# Units constants
kelvin = 1
picoseconds = 1
nanometer = 1

def printFlexLine(line, joint):
	return ("%d %d %s %s %s %s %d %s %s %s %d %s %.5f %.5f\n" % \
(int(line[0]), int(line[1]), joint, line[3], line[4], line[5], int(line[6]), \
line[7], line[8], line[9], int(line[10]), line[11], float(line[12]), \
float(line[13]))
	)
#

class Platform:
	"""
	Contains information on # of threads and openmm GPU usage
	"""
	def __init__(self):
		pass
	#
	
	def getPlatformByName(nameOfPlatform):
		print(nameOfPlatform)
		if nameOfPlatform == 'GPU':
			return True
		else:
			return False
	#			
#

class PDBReporter:
	"""
	Reporter
	"""
	def __init__(self, outputDir, stride):
		self.outputDir = outputDir
		self.stride = stride
	#		
#

class World:
	""" 
	Contains the complement of the constraints. 
	In Gibbs sampling terminology: it contains specifications of 
	the blocks that are sampled.
	"""
	def __init__(self):
		pass
	#		
#

class System:
	""" Contains Worlds and forces """
	def __init__(self):
		self.moldirs = []
		self.topologyFNs = []
	#
	
	def addConstraints(self, FN):
		""" Read the flexible bonds from a file """
		pass
	#	
#

class AmberPrmtopFile:
	""" Read an Amber prmtop file """
	def __init__(self, FN):
		self.topology = None
		self.topologyFNs = FN
	#

	def createSystem(self, createDirs = True,
		nonbondedMethod = "NoCutoff",
        	nonbondedCutoff = 1.0*nanometer,
        	constraints = None,
        	rigidWater = True,
        	implicitSolvent = True,
        	soluteDielectric = 1.0,
        	solventDielectric = 78.5,
        	removeCMMotion = False
	):
		""" Create a Robosample system """
		system = System()

		if createDirs == True:
			# Create directory to store molecules
			moldirs = "robots"
			if not os.path.exists(moldirs):
				try:
					os.makedirs(moldirs)
				except OSError as error:
					if error.errno != errno.EEXIST:
						raise
	
			# Create directories for each molecule
			for i in range(99999):
				if not os.path.exists("robots/bot" + str(i)):
					moldirs += "/bot" + str(i)
					break

			# Check if the number of directories is too huge
			if moldirs == "robots":
				print("Error: there are already 99999 directories in robots/. Exiting...")
				exit(1)

			# 
			if not os.path.exists(moldirs):
				try:
					os.makedirs(moldirs)
				except OSError as error:
					if error.errno != errno.EEXIST:
						raise

			# Copy prmtop file into its bot directory
			shutil.copy2(self.topologyFNs, moldirs + str("/bot.prmtop"))

			# Create directory to store trajectories
			moldirs = "robots/pdbs"
			if not os.path.exists(moldirs):
				try:
					os.makedirs(moldirs)
				except OSError as error:
					if error.errno != errno.EEXIST:
						raise
	
		# Get topology files
		#topologyFNs = []
		#for dire in glob.glob("robots/bot*"):
		#	topologyFNs.append(os.listdir(dire))
		topologyFNs = glob.glob("robots/bot*/*.prmtop")

		#system.topologyFNs = list(itertools.chain.from_iterable(topologyFNs))
		for f in topologyFNs:
			path, filename = os.path.split(topologyFNs[0])
			system.moldirs.append(path)
			system.topologyFNs.append(filename)



		return system
	#
#

class AmberInpcrdFile:
	def __init__(self, FN):
		self.positionsFNs = FN

		self.positionsFs = open(self.positionsFNs, "r")
		all_lines = self.positionsFs.readlines()

		self.lines = []
		self.xyz = []
		j = -1
		k = -1
		for i in range(len(all_lines)):
			line = all_lines[i].rstrip('\n')
			j = j + 1
			if j >= 2:
				k = k + 1
				self.xyz.append([])
				self.xyz[k].append(float(line[ 0:12]))
				self.xyz[k].append(float(line[12:24]))
				self.xyz[k].append(float(line[24:36]))
				if len(line) >= 72:
					k = k + 1
					self.xyz.append([])
					self.xyz[k].append(float(line[36:48]))
					self.xyz[k].append(float(line[48:60]))
					self.xyz[k].append(float(line[60:72]))
					

		self.positions = np.array(self.xyz)
	#		
#

class Context:
	""" Current state """
	def __init__(self):
		self.positions = None
		self.positionsFNs = None
		self.mdtrajTraj = None
		self.path = None
		self.coils = []
		self.dihIxs = []
		self.sasa = None
		self.bondsReformat = []
		self.bondsBB = []
	#
	
	def setPositions(self, positions):
		self.positions = positions

		# Get last directory
		topologyFNs = glob.glob("robots/bot*/*.prmtop")
		if len(topologyFNs) == 0:
			print("No directory tree generated. Please use \"createDirs = True\" when creating system.")
			exit(1)
		if len(topologyFNs) > 99999:
			print("Too many topologies. Exiting...")
			exit(2)



		# Get all flexibility
		#lines = np.loadtxt(args.top)
		infile = open(topologyFNs[0], 'r')
		lines = infile.readlines()
		infile.close()
		
		flag = False
		dihlix = 0
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
						self.dihIxs.append(onedihIxs)
		self.dihIxs = np.array(self.dihIxs)
		

		# Get the last directory
		filename = None
		for f in topologyFNs:
			self.path, filename = os.path.split(f)
		
		inpTxt = """TITLE\n"""
		inpTxt += str(self.positions.shape[0]) + '\n'
		i = -1
		for pos in positions:
			i = i + 1
			if not i%2: # 0, 2, 4, ...
				inpTxt += ("%12.7f%12.7f%12.7f" % (pos[0], pos[1], pos[2]))
			else: # 1, 3, 5
				inpTxt += ("%12.7f%12.7f%12.7f\n" % (pos[0], pos[1], pos[2]))
		
		#inpTxt += '\n'

		self.positionsFNs = self.path + '/bot.rst7'
		inpF = open(self.positionsFNs, "w+")
		inpF.write(inpTxt)
		inpF.close()
	
		# Get molecule
		print(self.positionsFNs, topologyFNs)
		self.mdtrajTraj = md.load(self.positionsFNs, top = topologyFNs[0])
		
		# Get topology
		topology = self.mdtrajTraj.topology
		table, bonds = topology.to_dataframe()
		
		# Compute accesibility: shape (n_frames,n_atoms)
		self.sasa = md.shrake_rupley(self.mdtrajTraj, probe_radius = 0.1, mode = 'residue') # nm

		self.allflexFN = os.path.join(self.path, "bot.all.flex") # duplicate in simulation

		if not os.path.exists(self.allflexFN):
			allflexF = open(self.allflexFN, 'w')
			#print("#i1 i2 Joint # name1 elem1 resid1 resname1 name2 elem2 resid2 resname2 sasa1 sasa2")
			for i in range(topology.n_bonds):
				allflexF.write("%5d %5d Cartesian # %4s %s %6d %3s   %4s %s %6d %3s %8.5f %8.5f\n" % (bonds[i][0], bonds[i][1] \
					,table.values[ int(bonds[i][0]) ][1], table.values[ int(bonds[i][0]) ][2] \
					,table.values[ int(bonds[i][0]) ][3], table.values[ int(bonds[i][0]) ][4] \
					,table.values[ int(bonds[i][1]) ][1], table.values[ int(bonds[i][1]) ][2] \
					,table.values[ int(bonds[i][1]) ][3], table.values[ int(bonds[i][1]) ][4] \
					,self.sasa[0][ table.values[ int(bonds[i][0]) ][3] ], self.sasa[0][ table.values[ int(bonds[i][1]) ][3] ]
				))
			allflexF.close()
		###

		# All processing is done with MDTraj for now
		#self.mdtrajTraj = md.load(self.positionsFNs, top = topologyFNs[0])
		print("Calculating distance matrix...")
		self.distMat = md.compute_contacts(self.mdtrajTraj, contacts='all', scheme='ca', ignore_nonprotein=False, periodic=False)
		print("Calculating dssp...")
		self.dssp = md.compute_dssp(self.mdtrajTraj)
		print('dssp', self.dssp)

		self.coils = []
		coilsLen = 0
		currResIx = 0
		nextResIx = 1
		while nextResIx < len(self.dssp[0]):
			# Coil found
			if (self.dssp[0][currResIx] != 'C') and (self.dssp[0][nextResIx] == 'C'):
				self.coils.append([nextResIx, nextResIx + 1])
				coilsLen += 1

				#print(''.join(self.dssp[0][0:currResIx]), self.dssp[0][currResIx], ''.join(self.dssp[0][currResIx+1:]))
				#print('hc', self.coils)
			# Inside a coil
			elif (self.dssp[0][currResIx] == 'C') and (self.dssp[0][nextResIx] == 'C'):
				if (len(self.coils) == 0):
					self.coils.append([0, 1])
					coilsLen += 1
				else:
					self.coils[coilsLen - 1][1] += 1

				#print(''.join(self.dssp[0][0:currResIx]), self.dssp[0][currResIx], ''.join(self.dssp[0][currResIx+1:]))
				#print('cc', self.coils)
			# Coil end
			elif (self.dssp[0][currResIx] == 'C') and (self.dssp[0][nextResIx] != 'C'):
				#print(''.join(self.dssp[0][0:currResIx]), self.dssp[0][currResIx], ''.join(self.dssp[0][currResIx+1:]))
				#print('ch', self.coils)
				pass
			else:
				#print(''.join(self.dssp[0][0:currResIx]), self.dssp[0][currResIx], ''.join(self.dssp[0][currResIx+1:]))
				#print('hh', self.coils)
				pass

			currResIx += 1
			nextResIx += 1
		self.coils = np.array(self.coils)
		print('coils', self.coils)

		# Compute dot products between bonds
		for i in range(self.mdtrajTraj.topology.n_bonds):
			self.bondsReformat.append([int(bonds[i][0]), int(bonds[i][1])
				, table.values[ int(bonds[i][0]) ][1], table.values[ int(bonds[i][0]) ][3] 
				, table.values[ int(bonds[i][1]) ][1], table.values[ int(bonds[i][1]) ][3]
				, np.array(positions[int(bonds[i][0])] - positions[int(bonds[i][1])]) ]
			)
		#print('bondsReformat', self.bondsReformat)
	
		# Extract backbone
		for bvi in range(len(self.bondsReformat)):
			if (((self.bondsReformat[bvi][2] == 'N') and (self.bondsReformat[bvi][4] == 'CA')) or
			    ((self.bondsReformat[bvi][2] == 'CA') and (self.bondsReformat[bvi][4] == 'N')) or
			    ((self.bondsReformat[bvi][2] == 'CA') and (self.bondsReformat[bvi][4] == 'C')) or
			    ((self.bondsReformat[bvi][2] == 'C') and (self.bondsReformat[bvi][4] == 'CA'))):
				self.bondsBB.append(self.bondsReformat[bvi])
		
		# Actual calculation
		bondDots = np.empty((len(self.bondsBB), len(self.bondsBB))) * np.nan
		for bvi in range(len(self.bondsBB)):
			#print(self.bondsBB[bvi], end = ' ')
			for bvj in range(len(self.bondsBB)):
				#if bvj < bvi:
				if True:
					#print(bvj, end = ' ')
					v1 = self.bondsBB[bvi][6]
					v2 = self.bondsBB[bvj][6]
					n1 = np.linalg.norm(v1)
					n2 = np.linalg.norm(v2)
					bondDots[bvi, bvj] = bondDots[bvj, bvi] = np.abs(np.dot(v1, v2) / (n1*n2))
					print(bondDots[bvi, bvj], end = ' ')
			print()
		
		# Perform clustering on dot products
		import scipy.cluster.hierarchy as sch
		from scipy.spatial.distance import squareform
		from sklearn.cluster import AgglomerativeClustering

		reduced_distances = squareform(bondDots, checks=False) # Memory error
		Z = sch.linkage(reduced_distances, method='ward') # 'single', 'complete', 'weighted', and 'average'
		flatClusts = sch.fcluster(Z, 0.999, criterion='distance')
		flatClusts = np.array(flatClusts)
		print('flatClusts', flatClusts)
		#for clustNo in range(np.max(flatClusts)):
		for clustNo in range(1,2):
			for doti in range(flatClusts.shape[0]):
				if flatClusts[doti] == clustNo:
					#print(self.bondsBB[doti][0], self.bondsBB[doti][1], "Pin")
					print(doti, self.bondsBB[doti], self.bondsBB[doti])
		
		

		####

		#


	def process_flex(self, subset='rama', residRange=[0, 1], accRange=[0.0, 10.0], jointType='Pin', worldNo=2, FN=None):
		# Robosample tools
		from ls_parsetxt import ParseTxt
		
		# Read data
		pars = ParseTxt()
		pars.Read(self.allflexFN)
		pdata = pars.parsed_data
		
		aa = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
		sugars = ['UYB', '4YB', 'VMB', '0MA', '2MA', '0LB', '3LB', '6LB', '0fA', '0SA', 'YMA']

		if FN == None:
			FN = os.path.join(self.path, "bot." + str(worldNo) + ".flex")
		f = open(FN, "a+")

		# Get all
		if (subset).lower() == "all":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						f.write(printFlexLine(line, jointType))
		
		# Get phi
		if (subset).lower() == "phi":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						if(line[7] != "PRO"):
							if(((line[4] == "N") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "N"))):
								f.write(printFlexLine(line, jointType))
		
		# Get psi
		if (subset).lower() == "psi":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						if(((line[4] == "C") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "C"))):
							f.write(printFlexLine(line, jointType))
		
		# Get Ramachandran flexibility
		if (subset).lower() == "rama":
		#if "rama" in [sub.lower() for sub in subset]:
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						if(line[7] != "PRO"):
							if(((line[4] == "N") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "N"))):
								f.write(printFlexLine(line, jointType))
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						if(((line[4] == "C") and (line[8] == "CA")) or ((line[4] == "CA") and (line[8] == "C"))):
							f.write(printFlexLine(line, jointType))
		
		
		# Get sidechains
		if (subset).lower() == "side":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					line = pdata[pi]
					atom1Acc = float(line[-2])
					atom2Acc = float(line[-1])
					if( ((atom1Acc >= accRange[0]) and (atom1Acc <= accRange[1])) or 
					    ((atom2Acc >= accRange[0]) and (atom2Acc <= accRange[1])) ):
						if(line[7] == "ALA"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "VAL"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG1")) or ((line[4] == "CG1") and (line[8] == "CB")) or 
							   ((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "LEU"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "CD1")) or ((line[4] == "CD1") and (line[8] == "CG")) or 
							   ((line[4] == "CG") and (line[8] == "CD2")) or ((line[4] == "CD2") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "ILE"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG1")) or ((line[4] == "CG1") and (line[8] == "CB")) or 
							   ((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG1") and (line[8] == "CD1")) or ((line[4] == "CD1") and (line[8] == "CG1"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "MET"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
		
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "CG") and (line[8] == "SD")) or ((line[4] == "SD") and (line[8] == "CG"))):
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "CG") and (line[8] == "SD")) or ((line[4] == "SD") and (line[8] == "CG"))):
									f.write(printFlexLine(line, jointType))
		
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "SD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "SD"))):
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "SD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "SD"))):
									f.write(printFlexLine(line, jointType))
		
						if(line[7] == "PHE"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
						if(line[7] == "TRP"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
						if(line[7] == "TYR"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
		
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "CZ") and (line[8] == "OH")) or ((line[4] == "OH") and (line[8] == "CZ"))): 
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "CZ") and (line[8] == "OH")) or ((line[4] == "OH") and (line[8] == "CZ"))): 
									f.write(printFlexLine(line, jointType))
		
						if(line[7] == "ASN"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
						if(line[7] == "SER"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
		
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "CB") and (line[8] == "OG")) or ((line[4] == "OG") and (line[8] == "CB"))): 
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "CB") and (line[8] == "OG")) or ((line[4] == "OG") and (line[8] == "CB"))): 
									f.write(printFlexLine(line, jointType))
						if(line[7] == "THR"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG2")) or ((line[4] == "CG2") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "CB") and (line[8] == "OG1")) or ((line[4] == "OG1") and (line[8] == "CB"))):
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "CB") and (line[8] == "OG1")) or ((line[4] == "OG1") and (line[8] == "CB"))):
									f.write(printFlexLine(line, jointType))
						if(line[7] == "CYS"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if jointType in ["BallM", "BallF"] :
								if(((line[4] == "CB") and (line[8] == "SG")) or ((line[4] == "SG") and (line[8] == "CB"))): 
									f.write(printFlexLine(line, "Pin"))
							else:
								if(((line[4] == "CB") and (line[8] == "SG")) or ((line[4] == "SG") and (line[8] == "CB"))): 
									f.write(printFlexLine(line, jointType))
						if(line[7] == "GLN"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "LYS"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CD") and (line[8] == "CE")) or ((line[4] == "CE") and (line[8] == "CD"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CE") and (line[8] == "NZ")) or ((line[4] == "NZ") and (line[8] == "CE"))):
								f.write(printFlexLine(line, jointType))
						if(line[7] == "ARG"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CD") and (line[8] == "NE")) or ((line[4] == "NE") and (line[8] == "CD"))):
								f.write(printFlexLine(line, jointType))
							#if(((line[4] == "NE") and (line[8] == "CZ")) or ((line[4] == "CZ") and (line[8] == "NE"))):
							#	f.write(printFlexLine(line, jointType))
							#if(((line[4] == "CZ") and (line[8] == "NH1")) or ((line[4] == "NH1") and (line[8] == "CZ"))):
							#	f.write(printFlexLine(line, jointType))
							#if(((line[4] == "CZ") and (line[8] == "NH2")) or ((line[4] == "NH2") and (line[8] == "CZ"))):
							#	f.write(printFlexLine(line, jointType))
						if((line[7][0] == "H") and (line[7][1] == "I")):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
						if(line[7] == "ASP"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
						if(line[7] == "GLU"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "CD")) or ((line[4] == "CD") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
								
			
		# Glycan part	
		if (subset).lower() == "sugnln":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					if(line[7] == "NLN") and (line[11] == "NLN"):
							if(((line[4] == "CA") and (line[8] == "CB")) or ((line[4] == "CB") and (line[8] == "CA"))):
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CB") and (line[8] == "CG")) or ((line[4] == "CG") and (line[8] == "CB"))): 
								f.write(printFlexLine(line, jointType))
							if(((line[4] == "CG") and (line[8] == "ND2")) or ((line[4] == "ND2") and (line[8] == "CG"))):
								f.write(printFlexLine(line, jointType))
					if(line[7] == "NLN") and (line[11] in sugars):
						if((line[4] == "ND2") or (line[8] == "ND2")):
							f.write(printFlexLine(line, jointType))
			
		if (subset).lower() == "suginter":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					if(line[7] in sugars) and (line[11] in sugars) and (line[7] != line[11]):
						f.write(printFlexLine(line, jointType))
		
		if (subset).lower() == "sugout":
			for pi in range(len(pdata)):
				line = pdata[pi]
				if((int(line[6]) >= residRange[0]) and (int(line[6]) <= residRange[1])):
					if(line[7] in sugars) and (line[7] == line[11]): # within the same sugar
						if(line[7] not in ['0SA']):
							if(((line[4] in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']) and (line[8] not in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5'])) or \
							   ((line[8] in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']) and (line[4] not in ['C1', 'C2', 'C3', 'C4', 'C5', 'O5']))):
								if((line[5] != 'H') and (line[9] != 'H')):
									f.write(printFlexLine(line, jointType))
						elif(line[7] in ['0SA']):
							if(((line[4] in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']) and (line[8] not in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6'])) or \
							   ((line[8] in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']) and (line[4] not in ['C2', 'C3', 'C4', 'C5', 'C6', 'O6']))):
								if((line[5] != 'H') and (line[9] != 'H')):
									f.write(printFlexLine(line, jointType))

		wroteToFile = f.tell()
		f.close()

		if wroteToFile == 0:
			return False
		else:
			return True
		
		
	class NMAtoFlex:
		def __init__(self, MDTrajObj, iModExec, tempRoot="./temp", Modes=1):
			
			"""
			Class that takes an MDTraj Trajectory object, ideally a prmtop/inpcrd pair 
			so the atom indices don't change, and runs iMod on the coordinates, then
			generates a flex file compatible with the Robosample package.
			
			Parameters
			----------
			
			MDTrajObj:  MDTrajTrajectory Object
			            Trajectory object for which the NMA will be done.
			iModExec:   str
			            Path to iMod executable (imode_gcc)
			tempRoot:   str
			            Path where temporary iMod-generated files will be stored

			Modes:		int
							How many modes to compute. Default = 20.

			"""
			
			## Slice any non-protein atoms from the MDTraj Trajectory Object
			ProteinAtoms = MDTrajObj.topology.select("protein")
			MDTrajObj = MDTrajObj.atom_slice(ProteinAtoms)
			
			## Generate a PDB
			PDBName = "{}.pdb".format(tempRoot)
			MDTrajObj.save_pdb(PDBName, force_overwrite=True)
		
			## Run the actual iMod command
			iModCommand = [iModExec, PDBName, "-o {}".format(tempRoot), "-x",
								"--save_fixfile", "-n {}".format(Modes)]
			iMod_sub = subprocess.Popen(iModCommand, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
			iMod_sub.communicate()
       
			## Generate Appropriate Arrays containing information
			## Parse the fixfile, which is the same regardless of mode
			flexfields = [("ResIx", int), ("Phi", float), ("Chi", float), ("Psi", float)]
			FixIn = open("{}.fix".format(tempRoot), "r").readlines()
			FixArray = np.zeros(0, dtype=flexfields)
			for Line in FixIn:
				LineSplit = Line.strip("\n").split()
				if (len(LineSplit) == 4):
					FixArray = np.append(FixArray,np.array([(LineSplit[0], LineSplit[1], LineSplit[2], LineSplit[3])], dtype=FixArray.dtype)) 

			## Parse the Evec file, generating a list of lists, one for each mode
			EvecIn = open("{}_ic.evec".format(tempRoot), "r").readlines()
			for lineIx in range(len(EvecIn)):
				EvecIn[lineIx] = EvecIn[lineIx].split()
        
			## Get the number of computed modes
			NOfModes = int(EvecIn[1][3])
			print ("Computed a number of {} modes".format(NOfModes))
			EvecIn = EvecIn[2:]
			EvecModes = []
			Frequencies = []
			
			## We retrieve the SIGNED values of each evec contribution
			CurrentMode=[]
			for lineIx in range(len(EvecIn)):
				if (EvecIn[lineIx][0] == "****"):
					if (len(CurrentMode) > 0):
						EvecModes.append(CurrentMode)
						CurrentMode=[]
					continue
				if (len(EvecIn[lineIx][0]) > 2):
					for EvecValue in EvecIn[lineIx]:
						CurrentMode.append(float(EvecValue)) 
				else:
					Frequencies.append(1/float(EvecIn[lineIx][1]))
			EvecModes.append(CurrentMode)

			EvecModesScale = np.array(EvecModes)

			self.EvecModes = EvecModes
			self.EvecModesScale = EvecModesScale
			self.MDTrajObj = MDTrajObj
			self.PDBName = PDBName
			self.FixArray = FixArray
			self.tempRoot = tempRoot

		def GetFlexScaled(self, Output="NMAFlexFiles/receptor", GenerateVMD=True):
			"""
			Function that uses the previously generated arrays to generate the
			Robosample-compatible flex files, along with a scaling factor for velocity.
			
			Parameters
			----------
			
			Output:     	str
			            Root name for generated flex files.
			
			GenerateVMD:	bool
							If True: generates a PDB file where the beta column is replaced by
							the NMA-based velocity scale 
							Default: True
			"""

			##Create Output dir if it doesn't exist:
			if (os.path.exists(Output.split("/")[0]) == False):
				print("{} did not exist and has been created.".format(Output))
				os.mkdir(Output.split("/")[0])

			ModeArray = self.FixArray.copy()
		
			##Assign vectors to dihedral angles and scaling factors
			for EvecModesScaleIx in range(len(self.EvecModes)):
				counter = 0
				for aIx in range(ModeArray.shape[0]):
					if (ModeArray[aIx]["Phi"] != 0):
						ModeArray[aIx]["Phi"] = self.EvecModesScale[EvecModesScaleIx][counter]
						counter +=1
					if (ModeArray[aIx]["Chi"] != 0):
						ModeArray[aIx]["Chi"] = self.EvecModesScale[EvecModesScaleIx][counter]
						counter +=1
					if (ModeArray[aIx]["Psi"] != 0):
						ModeArray[aIx]["Psi"] = self.EvecModesScale[EvecModesScaleIx][counter]
						counter +=1
	
				## Sort Arrays by 
				sortFlexChi = np.sort(ModeArray, order="Chi")
				sortFlexChi = np.flip(sortFlexChi)
				sortFlexPhi = np.sort(ModeArray, order="Phi")
				sortFlexPhi = np.flip(sortFlexPhi)
				sortFlexPsi = np.sort(ModeArray, order="Psi")
				sortFlexPsi = np.flip(sortFlexPsi)
			
				## Extract the associated PHI/PSI angles from the above array  
				fixfields = [("atom1", int), ("atom2", int),("scale", float)]
				joint="Pin" 
				FlexOutArray = np.zeros(0, dtype=fixfields)
				for ResIx in range(ModeArray.shape[0]):
					if (str(self.MDTrajObj.topology.residue(sortFlexPhi[ResIx]["ResIx"]))[0:3] != "PRO"):
						NIndex = self.MDTrajObj.topology.select("resid {} and name N".format(sortFlexPhi[ResIx]["ResIx"]))[0]
						CAIndex = self.MDTrajObj.topology.select("resid {} and name CA".format(sortFlexPhi[ResIx]["ResIx"]))[0]
						ScaleFactor = sortFlexPhi[ResIx]["Phi"]
						FlexOutArray = np.append(FlexOutArray,np.array((NIndex, CAIndex, ScaleFactor), dtype=fixfields)) 
				for ResIx in range(ModeArray.shape[0]):
					CAIndex = self.MDTrajObj.topology.select("resid {} and name CA".format(sortFlexPsi[ResIx]["ResIx"]))[0]
					CIndex = self.MDTrajObj.topology.select("resid {} and name C".format(sortFlexPsi[ResIx]["ResIx"]))[0]
					ScaleFactor = sortFlexPsi[ResIx]["Psi"]
					FlexOutArray = np.append(FlexOutArray,np.array((CAIndex, CIndex, ScaleFactor), dtype=fixfields)) 
				for ResIx in range(ModeArray.shape[0]):
					if (str(self.MDTrajObj.topology.residue(sortFlexChi[ResIx]["ResIx"]))[0:3] not in ["PRO", "GLY","ALA"]):
						CAIndex = self.MDTrajObj.topology.select("resid {} and name CA".format(sortFlexChi[ResIx]["ResIx"]))[0]
						CBIndex = self.MDTrajObj.topology.select("resid {} and name CB".format(sortFlexChi[ResIx]["ResIx"]))[0]
						ScaleFactor = sortFlexChi[ResIx]["Chi"]
						FlexOutArray = np.append(FlexOutArray,np.array((CAIndex, CBIndex, ScaleFactor), dtype=fixfields)) 
	
				## We now sort the above array by the ScaleFactor, and write a flex file
				FlexOutArray = np.sort(FlexOutArray, order="scale")
				FlexOutArray = np.flip(FlexOutArray)
				
				FlexOut = open("{}.Mode{}.Scaled.flex".format(Output, EvecModesScaleIx), "w")
				FlexOut.write("##Scaled Flex File##\n")
				for Ix in range(FlexOutArray.shape[0]):
					if (FlexOutArray[Ix]["scale"] != 0):
						FlexOut.write("{} {} {} {} #{} {} \n".format(FlexOutArray[Ix]["atom1"], FlexOutArray[Ix]["atom2"], joint, FlexOutArray[Ix]["scale"], self.MDTrajObj.topology.atom(FlexOutArray[Ix]["atom1"]), self.MDTrajObj.topology.atom(FlexOutArray[Ix]["atom2"])))
				FlexOut.close()

				## Write PDB with beta factors that reflect the scales
				if (GenerateVMD == True):
					BetaAtoms = list(FlexOutArray["atom1"])
					BetaFactors = list(FlexOutArray["scale"])
					DefaultBeta = -1 
	
					BetaColor = []
					BetaColor.append(BetaAtoms)
					BetaColor.append(BetaFactors)
	
					PDBFileIn = open(self.PDBName, "r")				
					PDBLines = []
					for line in PDBFileIn:
						PDBLines.append(line.strip().split())
					for line in PDBLines:
						if (line[0] == "ATOM"):
							if (int(line[1])-1 in BetaColor[0]):
								line[10] = BetaColor[1][BetaColor[0].index(int(line[1])-1)]*100
							else:
								line[10] = DefaultBeta		
					PDBFileIn = open(self.PDBName, "w")
					for i in PDBLines:
						if (len(i) == 12):
							PDBFileIn.write("{:6s}{:5d}  {:^4s}{:3s} {:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}\n".format(i[0], int(i[1]), i[2], i[3], i[4], int(i[5]), float(i[6]), float(i[7]), float(i[8]), float(i[9]), float(i[10]), i[11]))	##PDB Format
					PDBFileIn.close()
	
					print("VMD File written!")

		def CleanUp(self):
			"""
			Simple function that cleans up the files generated during the iMod process.
			Does NOT remove the pdb generated, since you need that for the VMD scripts to work. 
			"""
			for FileExt in [".fix", ".log", "_ic.evec", "_model.pdb"]:
				procList   = ["rm", "{}{}".format(self.tempRoot,FileExt)]
				removeTemp = subprocess.Popen(procList, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	#	
#

class Simulation:
	def __init__(self, topology, system, integrator='HMC', platform='CPU', properties={'nofThreads': 2}, addDefaultWorld=True):
		print("Starting Simulation init")
		self.context = Context()
		self.system = system
		self.integrator = integrator
		self.openmmTrue = platform

		self.nofThreads = properties['nofThreads']

		# Get last directory
		self.topologyFNs = glob.glob("robots/bot*/*.prmtop")
		if len(self.topologyFNs) > 99999:
			print("Too many topologies. Exiting...")
			exit(2)

		# Get the last directory
		path = None
		filename = None
		for f in self.topologyFNs:
			path, filename = os.path.split(f)
		self.path = path
		self.filename = filename
		
		print("DEBUG", self.system.topologyFNs)
		self.topFN = os.path.join(self.path, self.system.topologyFNs[0])
		self.crdFN = os.path.join(self.path, "bot.rst7")

		self.allflexFN = os.path.join(self.path, "bot.all.flex")
		self.defaultRB = os.path.join(self.path, "bot.rb")
		self.defaultFLEX = os.path.join(self.path, "bot.flex")


		# Get all bonds
		#if not os.path.exists(self.allflexFN):
			
			#execList = ['python3', '/home/pcuser/git3/Robosample/tools/getAllBonds.py', '--top', str(self.topFN), "--traj", self.crdFN, '--flex', 'bot.flex', '--probesize', '0.1']
			#proc = None
			#with open(os.path.join(self.path, "bot.all.flex"), "w") as outF:
			#	print("executing getAllBondsi with", self.crdFN)
			#	proc = subprocess.run(execList, env={**os.environ}, stdout = outF)
			
#		execStr = "python3 $ROBOSAMPLEDIR/tools/getAllBonds.py --top " + str(self.topFN) \
#			+ " --traj " + self.crdFN \
#			+ " --flex bot.flex --probesize 0.1 >" + os.path.join(self.path, "bot.all.flex")
#		os.system("echo \"" + execStr + "\"")
#		os.system(execStr)

		self.reporters = []
		self.inpDict = None
		self.nofWorlds = 0

		if addDefaultWorld == True:
			self.inpDict = {
				'MOLECULES': ['robots/bot0'],
				'ROUNDS': [2000],
				'DISTANCE': [0, 1],
				'DIHEDRAL': [0, 1, 2, 3],
				'OUTPUT_DIR': ['robots/'],
				'RANDOM_WORLD_ORDER': ['FALSE'], 
	
				'PRMTOP': ['bot.prmtop'],
				'INPCRD': ['bot.rst7'],
				'RBFILE': ['bot.rb'],
				'FLEXFILE': ['bot.all.flex'],
				'ROOT_MOBILITY': ['Cartesian'],
				'RUN_TYPE': ['Normal'],
				'ROUNDS_TILL_REBLOCK': [10],
				'WORLDS': ['R0'],
				'ROOTS': [0],
				'SAMPLER': [self.integrator.type],
				'TIMESTEPS': [0.0007],
				'MDSTEPS': [5],
				'BOOST_MDSTEPS': [1],
				'SAMPLES_PER_ROUND': [4],
				'REPRODUCIBLE': ['FALSE'],
				'SEED': [1],
				'THERMOSTAT': ['Andersen'],
				'TEMPERATURE_INI': [self.integrator.T],
				'TEMPERATURE_FIN': [self.integrator.T],
				'BOOST_TEMPERATURE': [600],
				'FFSCALE': ['AMBER'],
				'GBSA': [1.0],
				'FIXMAN_POTENTIAL': ['FALSE'],
				'FIXMAN_TORQUE': ['FALSE'],
				'VISUAL': ['FALSE'],
				'PRINT_FREQ': [1],
				'WRITEPDBS': [1],
				'GEOMETRY': ['FALSE'],
				'THREADS': [self.nofThreads],
				'OPENMM': [str(self.openmmTrue).upper()]
			}
	
			self.nofWorlds = 1
		else:
			self.inpDict = {
				'MOLECULES': ['robots/bot0'],
				'ROUNDS': [2000],
				'DISTANCE': [0, 1],
				'DIHEDRAL': [0, 1, 2, 3],
				'OUTPUT_DIR': ['robots/'],
				'RANDOM_WORLD_ORDER': ['FALSE'], 
	
				'PRMTOP': [],
				'INPCRD': [],
				'RBFILE': [],
				'FLEXFILE': [],
				'ROOT_MOBILITY': [],
				'RUN_TYPE': [],
				'ROUNDS_TILL_REBLOCK': [],
				'WORLDS': [],
				'ROOTS': [],
				'SAMPLER': [],
				'TIMESTEPS': [],
				'MDSTEPS': [],
				'BOOST_MDSTEPS': [],
				'SAMPLES_PER_ROUND': [],
				'REPRODUCIBLE': [],
				'SEED': [],
				'THERMOSTAT': [],
				'TEMPERATURE_INI': [],
				'TEMPERATURE_FIN': [],
				'BOOST_TEMPERATURE': [],
				'FFSCALE': [],
				'GBSA': [],
				'FIXMAN_POTENTIAL': [],
				'FIXMAN_TORQUE': [],
				'VISUAL': [],
				'PRINT_FREQ': [],
				'WRITEPDBS': [],
				'GEOMETRY': [],
				'THREADS': [],
				'OPENMM': []
			}
	
			self.nofWorlds = 0

		print("Done Simulation init")
	#

	def inpEverythingButFlex(self, rootMobility='Weld', timestep=0.001, mdsteps=10, contacts=False, contactCutoff=0, samples=1):
		self.inpDict['PRMTOP'].append('bot.prmtop')
		self.inpDict['INPCRD'].append('bot.rst7')
		self.inpDict['RBFILE'].append('bot.rb')
		self.inpDict['ROOT_MOBILITY'].append(rootMobility)
		self.inpDict['RUN_TYPE'].append('Normal')
		self.inpDict['ROUNDS_TILL_REBLOCK'].append(10)
		self.inpDict['ROOTS'].append(0)
		self.inpDict['SAMPLER'].append(self.integrator.type)
		self.inpDict['TIMESTEPS'].append(timestep)
		self.inpDict['MDSTEPS'].append(mdsteps)
		self.inpDict['BOOST_MDSTEPS'].append(1)
		self.inpDict['SAMPLES_PER_ROUND'].append(samples)
		self.inpDict['REPRODUCIBLE'].append('FALSE')
		self.inpDict['SEED'].append(1)
		self.inpDict['THERMOSTAT'].append('Andersen')
		self.inpDict['TEMPERATURE_INI'].append(self.integrator.T)
		self.inpDict['TEMPERATURE_FIN'].append(self.integrator.T)
		self.inpDict['BOOST_TEMPERATURE'].append(600)
		self.inpDict['FFSCALE'].append('AMBER')
		self.inpDict['GBSA'].append(1.0)
		self.inpDict['FIXMAN_POTENTIAL'].append('TRUE')
		self.inpDict['FIXMAN_TORQUE'].append('TRUE')
		self.inpDict['VISUAL'].append('FALSE')
		self.inpDict['PRINT_FREQ'].append(1)
		if self.inpDict['WRITEPDBS'] == []:
			self.inpDict['WRITEPDBS'].append(1)
		else:
			self.inpDict['WRITEPDBS'].append(0)
		self.inpDict['GEOMETRY'].append('FALSE')
		self.inpDict['THREADS'].append(self.nofThreads)
		self.inpDict['OPENMM'].append(str(self.openmmTrue).upper())

		self.inpDict['WORLDS'].append('R' + str(self.nofWorlds))

	def addWorld(self, regionType='stretch', region=[[1, 2]], rootMobility='Weld', 
		timestep=0.001, mdsteps=10, argJointType="Pin", subsets=['rama'], 
		contacts=False, contactCutoff=0,
		samples=1):

		print("Starting Simulation addWorld")
		self.nofWorlds += 1

		region = np.array(region)
		if region.ndim != 2:
			print("Error adding world. You need to specify a 2-dim array of stretches.")
			exit(3)

		nofAAs = 0
		for reg in region:
			nofAAs += (reg[1] - reg[0])

		if regionType == 'roll':
			# Recursive call stretch like
			for regi in range(region.shape[0]):
				for regj in range(region[regi][0], region[regi][1] + 1):
					# Add specified region
					outFN = os.path.join(self.path, "bot.bb." + str(regj) + ".phi.flex")
	
					if os.path.isfile(outFN):
						os.remove(outFN)
	
					wroteToFile = False
					for subsetSpec in subsets:
						process_flexReturn = self.context.process_flex(subset='phi',
							residRange = [regj, regj],
							accRange=[0.0, 10.0],
							jointType=argJointType,
							worldNo=0, FN = outFN)
						wroteToFile = wroteToFile or process_flexReturn
	
					# Add flex filename to Robosample input
					if wroteToFile == True:
						dirName, fileName = os.path.split(outFN)
						self.inpDict['FLEXFILE'].append(fileName)
						self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)
						print("Wrote to", fileName)
		
						self.nofWorlds += 1
	
					# Add specified region
					outFN = os.path.join(self.path, "bot.bb." + str(regj) + ".psi.flex")
	
					if os.path.isfile(outFN):
						os.remove(outFN)
	
					wroteToFile = False
					for subsetSpec in subsets:
						process_flexReturn = self.context.process_flex(subset='psi',
							residRange = [regj, regj],
							accRange=[0.0, 10.0],
							jointType=argJointType,
							worldNo=0, FN = outFN)
						wroteToFile = wroteToFile or process_flexReturn
	
					# Add flex filename to Robosample input
					if wroteToFile == True:
						dirName, fileName = os.path.split(outFN)
						self.inpDict['FLEXFILE'].append(fileName)
						self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)
						print("Wrote to", fileName)
		
						self.nofWorlds += 1

		#

		# 
		elif regionType == 'stretch':
			outFN = os.path.join(self.path, "bot.stretch." + str(self.nofWorlds) + ".flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			# Add specified region
			for i in range(0, region.shape[0]):
				for subsetSpec in subsets:
					self.context.process_flex(subset=subsetSpec,
						residRange=region[i],
						accRange=[0.0, 10.0],
						jointType=argJointType,
						worldNo=0, FN = outFN)

			# Add contacts if specified
			print("region", region)
			flatRegion = []
			for rangePair in region:
				flatRegion.append(range(rangePair[0], rangePair[1]+1))
			flatRegion = np.array(flatRegion).flatten()
			print("flatRegion", flatRegion)

			contactList = []
			pi = -1
			for pair in self.context.distMat[1]:
				pi += 1
				if (pair[0] in flatRegion) or (pair[1] in flatRegion):
					if self.context.distMat[0][0][pi] < contactCutoff:
						if pair[0] not in flatRegion:
							contactList.append([pair[0], pair[0]])
						if pair[1] not in flatRegion:
							contactList.append([pair[1], pair[1]])

			contactList = np.array(contactList)
			#print("addWorld contacts:", contactList)
			for i in range(0, contactList.shape[0]):
				self.context.process_flex(subset='side',
					residRange=contactList[i],
					accRange=[0.1, 10.0],
					jointType='Pin',
					worldNo=0, FN = outFN)


			# Add flex filename to Robosample input
			dirName, fileName = os.path.split(outFN)
			self.inpDict['FLEXFILE'].append(fileName)
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)
			print("Wrote to", fileName)
			
		# Get acc, loops and sugars
		elif regionType in ['coils', 'loops']:
			outFN = os.path.join(self.path, "bot.loops.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for i in range(0, self.context.coils.shape[0]):
				print("Adding to flex region", self.context.coils[i])
				for subsetSpec in subsets:
					self.context.process_flex(subset=subsetSpec,
						residRange=self.context.coils[i],
						accRange=[0.0, 10.0],
						jointType=argJointType,
						worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.loops.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		# Get acc, loops and sugars
		elif regionType == 'accesible':
			outFN = os.path.join(self.path, "bot.acc.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for subsetSpec in subsets:
				self.context.process_flex(subset=subsetSpec,
					residRange=[0, self.context.positions.shape[0]],
					accRange=[0.5, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

				self.context.process_flex(subset=subsetSpec,
					residRange=[0, self.context.positions.shape[0]],
					accRange=[0.5, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.acc.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		# Get acc, loops and sugars
		elif regionType == 'sugars':
			outFN = os.path.join(self.path, "bot.sugars.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='sugnln',
					residRange=[0, self.context.positions.shape[0]],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='suginter',
					residRange=[0, self.context.positions.shape[0]],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='sugout',
					residRange=[0, self.context.positions.shape[0]],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.sugars.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		# Get acc, loops and sugars
		elif regionType == 'sugnln':
			outFN = os.path.join(self.path, "bot.sugnln.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='sugnln',
					residRange=region[i],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.sugnln.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		elif regionType == 'suginter':
			outFN = os.path.join(self.path, "bot.suginter.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='suginter',
					residRange=region[i],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.suginter.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		elif regionType == 'sugout':
			outFN = os.path.join(self.path, "bot.sugout.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			for i in range(0, region.shape[0]):
				self.context.process_flex(subset='sugout',
					residRange=region[i],
					accRange=[0.0, 10.0],
					jointType=argJointType,
					worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.sugout.flex')
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)

		# Ball world		
		elif regionType == 'ball':
			outFN = os.path.join(self.path, "bot.ball.flex")
			if os.path.isfile(outFN):
				os.remove(outFN)

			self.context.process_flex(subset='rama',
				residRange=[0, self.context.positions.shape[0]],
				accRange=[0.0, 10.0],
				jointType=argJointType,
				worldNo=0, FN = outFN)

			self.context.process_flex(subset='side',
				residRange=[0, self.context.positions.shape[0]],
				accRange=[0.0, 10.0],
				jointType=argJointType,
				worldNo=0, FN = outFN)

			self.inpDict['FLEXFILE'].append('bot.ball.flex')
		
			self.inpEverythingButFlex(rootMobility, timestep, mdsteps, contacts, contactCutoff, samples)
		print("Done Simulation addWorld")
	#

	def step(self, nofSteps):
		print("Starting Simulation step")
		os.system("touch " + self.defaultRB)
		os.system("touch " + self.defaultFLEX)


		# Put input together
		inpTxt = '''# Robosample input \n'''

		moldirs = self.system.moldirs
		for moldir in self.system.moldirs:
			#inpTxt += (' ' + moldir)
			pass

		for key in list(self.inpDict.keys()):
			inpTxt += key
			for val in self.inpDict[key]:
				inpTxt +=  " " + str(val)
			inpTxt += '\n'
		
		# Modify input
		self.inpDict['ROUNDS'] = [nofSteps]
		self.inpDict['OUTPUT_DIR'] = [self.reporters[0].outputDir]

		# Write input file
		inpFN = 'inp.test'
		inpF = open(inpFN, "w+")
		inpF.write(inpTxt)
		inpF.close()

		#os.system("echo \'" + inpTxt + "\'")
		#os.system("$ROBOSAMPLEDIR/build?debug/tests/Robosample inp.test")

		print("Done Simulation step")
	#	
#


class HMCIntegrator:
	def __init__(self, T, ts):
		self.T = T
		self.ts = ts

		self.type = 'HMC'
	#		
#

class VVIntegrator:
	def __init__(self, T, ts):
		self.T = T
		self.ts = ts
		
		self.type = 'VV'
	#		
#

class LAHMCIntegrator:
	def __init__(self, T, ts):
		self.type = 'LAHMC'
	#		
#

class NUTSIntegrator:
	def __init__(self, T, ts):
		self.type = 'NUTS'
	#		
#


