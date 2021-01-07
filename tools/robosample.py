# Imports
import os, errno
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

	def createSystem(self, createDirs = True):
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
	def __init__(self, topology, system, integrator):
		self.context = Context()
		self.system = system
		self.integrator = integrator

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
	#

	def step(self, nofSteps):
		topFN = os.path.join(self.path, self.system.topologyFNs[0])
		crdFN = os.path.join(self.path, "bot.rst7")
		defaultRB = os.path.join(self.path, "bot.rb")
		defaultFLEX = os.path.join(self.path, "bot.flex")
		os.system("touch " + defaultRB)
		os.system("touch " + defaultFLEX)

		# All processing is done with MDTraj for now
		md.load(self.context.positionsFNs, top = topFN)

		# Get all bonds
		allflexFN = os.path.join(self.path, "bot.all.flex")
		if not os.path.exists(allflexFN):
			os.system("python3 $ROBOSAMPLEDIR/tools/getAllBonds.py --top " + str(topFN) 
				+ " --traj " + crdFN
				+ " --flex bot.flex --probesize 0.1 >" + os.path.join(self.path, "bot.all.flex")
			)

		# Get hinges, loops and sugars
		execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset rama --accRange 0.5 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
				+ " > " + os.path.join(self.path, "bot.hinges.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)
		execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset side --accRange 0.5 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
				+ " >> " + os.path.join(self.path, "bot.hinges.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)
		execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset sugnln --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
				+ " >> " + os.path.join(self.path, "bot.hinges.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)
		execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset suginter --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
				+ " >> " + os.path.join(self.path, "bot.hinges.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)
		execStr = "python3 $ROBOSAMPLEDIR/tools/process_flex.py --inFN " + os.path.join(self.path, "bot.all.flex") \
				+ " --subset sugout --accRange 0.0 10 --joint Pin --residRange 0 " + str(self.context.positions.shape[0]) \
				+ " >> " + os.path.join(self.path, "bot.hinges.flex")
		os.system("echo \"" + execStr + "\"")
		os.system(execStr)

		# Ball world		
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

		# Put input together
		inpTxt = '''# Robosample input \n'''

		moldirs = self.system.moldirs
		inpTxt += ("MOLECULES")
		for moldir in self.system.moldirs:
			inpTxt += (' ' + moldir)

		restOfTxt = ''' 
PRMTOP bot.prmtop bot.prmtop       # Parameter file
INPCRD bot.rst7 bot.rst7     # Coordinate / Restart file
RBFILE bot.rb bot.rb   # Rigid bodies definition file
FLEXFILE bot.all.flex bot.hinges.flex # Flexibility definition file
ROOT_MOBILITY Cartesian Free # Ground to Compound mobilizer
OUTPUT_DIR robots robots

# Simulation
RUN_TYPE Normal Normal # normal HMC or Non-Eq HMC
ROUNDS ''' + str(nofSteps) + '''
ROUNDS_TILL_REBLOCK 10 10
RANDOM_WORLD_ORDER FALSE FALSE
WORLDS R0 R1                        # Regimen (IC, TD, MIX, RB, RBMIX)
ROOTS 0 0
SAMPLER ''' + self.integrator.type + " " + self.integrator.type + '''
TIMESTEPS 0.001 0.009               # Timesteps to be used with regimens
MDSTEPS 30 30  # Number of MD trial steps
BOOST_MDSTEPS 1 1
SAMPLES_PER_ROUND 1 3  # Number of MC trials within a mixing round
REPRODUCIBLE FALSE FALSE 
SEED 999 999

# Thermodynamics
THERMOSTAT Andersen Andersen   # Thermostat (HMC + Andersen = MD)
TEMPERATURE_INI  ''' + str(self.integrator.T) + " " + str(self.integrator.T) + '''
TEMPERATURE_FIN  ''' + str(self.integrator.T) + " " + str(self.integrator.T) + '''
BOOST_TEMPERATURE  1 1      # Temperature for constant temperature simulations
FFSCALE AMBER AMBER     # Force field scale factors
GBSA 0 0         # GBSA scale factor

# Generalized coordinates related
FIXMAN_POTENTIAL TRUE TRUE # Use Fixman potential
FIXMAN_TORQUE TRUE TRUE   # Use Fixman torque

# Output
VISUAL TRUE TRUE             # Use the visualizer
VISUAL FALSE FALSE             # Use the visualizer
PRINT_FREQ  1 1
WRITEPDBS 1 0    # Write pdbs
GEOMETRY FALSE FALSE         # Calculate geometric features

DISTANCE 1 2
DIHEDRAL 1 2 3 4 1 2 3 4

# Software specs
THREADS 0 0
OPENMM TRUE TRUE
''' 
		inpTxt += restOfTxt

		inpFN = 'inp.test'
		inpF = open(inpFN, "w+")
		inpF.write(inpTxt)
		inpF.close()

		#os.system("echo \'" + inpTxt + "\'")

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


