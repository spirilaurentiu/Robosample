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
	""" Conteins Worlds and forces """
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

			# Check if the number of directories is too hugh
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

		# Get the last directory
		path = None
		filename = None
		for f in topologyFNs:
			path, filename = os.path.split(f)
		
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


		self.positionsFNs = path + '/bot.rst7'
		inpF = open(self.positionsFNs, "w+")
		inpF.write(inpTxt)
		inpF.close()
	#	
#

class Simulation:
	def __init__(self, topology, system, integrator='HMC', platform='CPU', properties={'nofThreads': 2}):
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
		
		self.topFN = os.path.join(self.path, self.system.topologyFNs[0])
		self.crdFN = os.path.join(self.path, "bot.rst7")

		# All processing is done with MDTraj for now
		#md.load(self.context.positionsFNs, top = self.topFN)

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
			
		execStr = "python3 $ROBOSAMPLEDIR/tools/getAllBonds.py --top " + str(self.topFN) \
			+ " --traj " + self.crdFN \
			+ " --flex bot.flex --probesize 0.1 >" + os.path.join(self.path, "bot.all.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)

		self.reporters = []
		self.inpDict = None

		self.inpDict = {
			'MOLECULES': ['robots/bot0'],
			'ROUNDS': [10],
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
			'TIMESTEPS': [0.001],
			'MDSTEPS': [10],
			'BOOST_MDSTEPS': [1],
			'SAMPLES_PER_ROUND': [3],
			'REPRODUCIBLE': ['FALSE'],
			'SEED': [999],
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

		print("Done Simulation init")
	#

	def addWorld(self, regionType='stretch', region=[[1, 2]], rootMobility='Weld', timestep=0.001, mdsteps=10):
		print("Starting Simulation addWorld")
		self.nofWorlds += 1
		if regionType == 'stretch':
			region = np.array(region)
			if region.ndim != 2:
				print("Error adding world. You need to specify a 2-dim array of stretches.")
				exit(3)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset rama --accRange 0.0 10 --joint Pin --residRange " +  str(region[0][0]) + " "  + str(region[0][1]) \
				+ " > " + os.path.join(self.path, "bot.stretch" + str(self.nofWorlds) + ".flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset side --accRange 0.0 10 --joint Pin --residRange " +  str(region[0][0]) + " "  + str(region[0][1]) \
				+ " >> " + os.path.join(self.path, "bot.stretch" + str(self.nofWorlds) + ".flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			for i in range(1, region.shape[0]):
				execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset rama --accRange 0.0 10 --joint Pin --residRange " +  str(region[i][0]) + " "  + str(region[i][1]) \
					+ " >> " + os.path.join(self.path, "bot.stretch" + str(self.nofWorlds) + ".flex")
				os.system("echo \"" + execStr + "\"")
				os.system(execStr)
				execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset side --accRange 0.0 10 --joint Pin --residRange " +  str(region[i][0]) + " "  + str(region[i][1]) \
					+ " >> " + os.path.join(self.path, "bot.stretch" + str(self.nofWorlds) + ".flex")
				os.system("echo \"" + execStr + "\"")
				os.system(execStr)
			
			self.inpDict['FLEXFILE'].append("bot.stretch" + str(self.nofWorlds) + ".flex")
			
		# Get acc, loops and sugars
		if regionType == 'accesible':
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset rama --accRange 0.5 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " > " + os.path.join(self.path, "bot.acc.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset side --accRange 0.5 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " >> " + os.path.join(self.path, "bot.acc.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset sugnln --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " >> " + os.path.join(self.path, "bot.acc.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset suginter --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " >> " + os.path.join(self.path, "bot.acc.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset sugout --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " >> " + os.path.join(self.path, "bot.acc.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
	
			self.inpDict['FLEXFILE'].append('bot.acc.flex')

		# Ball world		
		if regionType == 'ball':
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset rama --accRange 0.0 10 --joint BallM --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " > " + os.path.join(self.path, "bot.ball.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)
			execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
					+ " --subset side --accRange 0.0 10 --joint BallM --residRange 0 " + str(self.context.positions.shape[0]) \
					+ " >> " + os.path.join(self.path, "bot.ball.flex")
			os.system("echo \"" + execStr + "\"")
			os.system(execStr)

			self.inpDict['FLEXFILE'].append('bot.ball.flex')

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
		self.inpDict['SAMPLES_PER_ROUND'].append(1)
		self.inpDict['REPRODUCIBLE'].append('FALSE')
		self.inpDict['SEED'].append(999)
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
		self.inpDict['WRITEPDBS'].append(0)
		self.inpDict['GEOMETRY'].append('FALSE')
		self.inpDict['THREADS'].append(self.nofThreads)
		self.inpDict['OPENMM'].append(str(self.openmmTrue).upper())

		self.inpDict['WORLDS'].append('R' + str(self.nofWorlds))
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

		os.system("echo \'" + inpTxt + "\'")
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


