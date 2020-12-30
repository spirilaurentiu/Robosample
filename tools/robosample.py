# Imports
import os

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
		pass
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
	#

	def createSystem(self):
		""" Create a Robosample system """
		return System()
	#
#

class AmberInpcrdFile:
	def __init__(self, FN):
		self.positions = None
	#		
#

class Context:
	""" Current state """
	def __init__(self):
		pass
	#
	
	def setPositions(self, positions):
		pass
	#	
#

class Simulation:
	def __init__(self, topology, system, integrator):
		self.context = Context()
	#

	def step(self, nofSteps):
		os.system('echo Stepping to ' + str(nofSteps))
	#	
#


class HMCIntegrator:
	def __init__(self, T, ts):
		pass
	#		
#

class VVIntegrator:
	def __init__(self, T, ts):
		pass
	#		
#

class LAHMCIntegrator:
	def __init__(self, T, ts):
		pass
	#		
#

class NUTSIntegrator:
	def __init__(self, T, ts):
		pass
	#		
#


