import sys, os, glob
import numpy as np
import numpy.linalg
import scipy 
import scipy.stats
from pymbar import timeseries as ts
from pymbar.utils import ParameterError


# Class that integrates several methods to compute autocorrelation data
class Autocorrelation:
	def __init__(self):
		# Input filename
		inFN = None
		stride = 1 # jump every stride line
		# nof columns to be analyzed
		ncols = None
		# Actual data
		cols = []
		tiny = 0.0000001
	#

	# Loads data from an input file
	def loadDataFromFile(self, argInFN, argSkip_header, argSkip_footer, argUsecols, argStride):
		inFN = argInFN
		self.stride = argStride
		self.ncols = len(argUsecols)

		with open(inFN, 'r') as in_F:
			#print "Input file", inFN
			alldata = np.genfromtxt(in_F, skip_header = argSkip_header, \
				skip_footer = argSkip_footer, usecols = argUsecols, \
				invalid_raise = False )
			alldata = alldata[::self.stride]
			print ("alldata.shape ", str(alldata.shape))
			print ("alldata ", alldata)
		
			# Print general stats
			print ("alldata meand stdev skew kurt ", np.mean(alldata), np.std(alldata), scipy.stats.skew(alldata) , scipy.stats.kurtosis(alldata))
		
			# Reshape data (transpose)
			if np.size(alldata.shape) == 1:
				self.cols = np.zeros((1, alldata.shape[0]))
				self.cols[0] = alldata
			else:
				shape_list = list(alldata.transpose().shape)
				shape_list[0] = shape_list[0] + 1
				shape_tuple = tuple(shape_list)
				self.cols = np.zeros(shape_tuple)
				for coli in range(self.ncols):
					self.cols[coli] = alldata.transpose()[coli]
		
			# Rehaped data
			#print "Reshaped data:"
			#for coli in range(self.ncols):
			#    print 'col', coli, 'shape', self.cols[coli].shape
			#    print 'col', coli, self.cols[coli]
	#

	def getData(self, argData):
		self.cols = argData
		self.ncols = argData.shape[0]
	#

	# Use NumPy correlate function
	def numpyCorr(self):
		results = np.zeros((self.ncols, self.cols[0].size))
		for coli in range(self.ncols):
			result = np.correlate(self.cols[coli], self.cols[coli], mode='full')
			results[coli] = result[result.size/2:]
		return results
	#

	# FFT
	def fftCorr(self):
		results = np.zeros((self.ncols, self.cols[0].size ))
		for coli in range(self.ncols):
			result = np.fft.ifft(np.abs(np.fft.fft(self.cols[coli]))**2).real
			results[coli] = result[:len(self.cols[coli])]
		return results

	# PyMBAR autocorrelation time
	def pymbarDetectEquilibration(self):
		results = np.zeros((self.ncols, 3))
		for coli in range(self.ncols):
			result = ts.detectEquilibration(self.cols[coli])
			results[coli] = result
		return results



	# Modified version of PyMBAR detectEquilibration
	# to be based on fft
	def pymbarDetectEquilibration_fft(self, fast=True, nskip=1):
		"""
		More documentation can be found in pymbar timeseries detectEquilibration function
		"""
		#results = np.zeros((self.ncols, 3))
		results = self.ncols * [4 * [None]]

		for coli in range(self.ncols):
			T = (self.cols[coli]).size
	   
			# Special case if timeseries is constant.
			if (self.cols[coli]).std() == 0.0:
				results[coli] = (0, 1, 1, np.ones([T - 1], np.float32))  # Changed from Neff=N to Neff=1 after issue #122
		
			g_t = np.ones([T - 1], np.float32)
			Neff_t = np.ones([T - 1], np.float32)
			for t in range(0, T - 1, nskip):
				try:
					g_t[t] = ts.statisticalInefficiency_fft(self.cols[coli][t:T])
				except ParameterError:  # Fix for issue https://github.com/choderalab/pymbar/issues/122
					g_t[t] = (T - t + 1)
				Neff_t[t] = (T - t + 1) / g_t[t]

			Neff_max = Neff_t.max()
			t = Neff_t.argmax()
			g = g_t[t]

			results[coli] = (t, g, Neff_max, Neff_t)
	
		return results
  
  
  
  
  
  
  
  
  
  
  
