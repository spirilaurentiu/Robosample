# Robosample tools
import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
from autocorFuncs import *
import Autocorrelation


class LogAnalyzer:
	def __init__(self, FNSeeds, simDirs, datacols, skip_header, skip_footer, stride):
		'''
		Initialize
		'''
		# Input parameters variables
		self.simDirs = simDirs
		self.FNSeeds = FNSeeds
		self.skip_header = skip_header
		self.skip_footer = skip_footer
		self.stride = stride
		self.nofDataCols = len(datacols)
		self.dataCols = datacols

		# Data storing variables
		self.logData = [None] * len(FNSeeds)
		self.trimLogData = [None] * len(FNSeeds)
		self.means = [None] * len(FNSeeds)
		self.variances = [None] * len(FNSeeds)
		self.stds = [None] * len(FNSeeds)
		self.mvAvgs = [None] * len(FNSeeds)
		self.avgAcc = [None] * len(FNSeeds)
		self.eqPoints = [None] * len(FNSeeds)
		self.eqPoint = [None] * len(FNSeeds)

		# Variables that store autocorrelation data
		self.fullAutocorrFunc = [None] * len(FNSeeds)
		self.trimAutocorrFunc = [None] * len(FNSeeds)
		self.cuts = [None] * len(FNSeeds)
		self.Iacs = [None] * len(FNSeeds)
		self.ESSs = [None] * len(FNSeeds)

		# PyMBAR based autocorrelation
		self.t_g_Neff = [None] * len(FNSeeds)
	#

	def Read(self, verbose = True):
		''' 
		Reads data from files and computes moments of series
		'''
		# Main loop: iterates through seeds = simulation reps
		for seedi in range(len(self.FNSeeds)):
			# Get file names
			FNList = []
			for di in range(len(self.simDirs)):
				FNList.append(glob.glob(os.path.join(self.simDirs[di] + "trim.log." + self.FNSeeds[seedi]))[0])
			nfiles = len(FNList)
		
			# Get raw data
			if verbose:
				print ("Gathering data from", FNList)
			rawdataChunks = []	
			for li in range(nfiles):
				with open(FNList[li], 'r') as in_FN1:
					rawdataChunks.append(np.genfromtxt(in_FN1, \
						skip_header=self.skip_header, skip_footer=self.skip_footer, \
						invalid_raise = False))
					rawdataChunks[li] = rawdataChunks[li][::self.stride]
			rawdataChunks = np.array(rawdataChunks)
			rawdata = np.concatenate(rawdataChunks)
		
			# Transpose data and add another column for acceptance
			self.logData[seedi] = np.concatenate((rawdata.transpose(), np.zeros((1, rawdata.shape[0]))), axis=0)
			ncols = self.logData[seedi].shape[0]
		
			# Attach acceptance indicator on the last column
			for i in range(rawdata.shape[0]):
				if rawdata[i][2] != rawdata[i][3]:
					self.logData[seedi][-1][i] = 1
			self.avgAcc[seedi] = np.mean(self.logData[seedi][-1])
	
			# Save moments for entire chosen timeseries
			self.variances[seedi] = np.zeros((self.nofDataCols))
			self.stds[seedi] = np.zeros((self.nofDataCols))
			self.means[seedi] = np.zeros((self.nofDataCols))
			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1
				self.means[seedi][dataColi] = np.mean(self.logData[seedi][coli])
				self.variances[seedi][dataColi] = np.var(self.logData[seedi][coli])
				self.stds[seedi][dataColi] = np.std(self.logData[seedi][coli])
		#
	
	def FindEquilibrationPoints(self, lagMaxima):
		'''
		Find the equilibration point for every data series
		'''
		for seedi in range(len(self.FNSeeds)):
			# Autocorrelation functions for data columns
			self.fullAutocorrFunc[seedi] = np.empty((self.nofDataCols, (self.logData[seedi][0]).size))
			maxAcor = 0
			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1
				# Autocorrelation time estimate (2009 Grossfield)
				M = lagMaxima[dataColi] # Calculate up to point M or until it reaches 0
				N = self.logData[seedi][coli].size # Sample size
				if M == 0: M = N
		
				# Autocorrelation functions calculated with different methods
				#self.fullAutocorrFunc[seedi][dataColi], cut, IAcFull, ESSfull = CestGrossfield(M, self.logData[seedi][coli])
				self.fullAutocorrFunc[seedi][dataColi], cut, IAcFull, ESSfull = autocorr1(M, self.logData[seedi][coli]) 
				#print "Full autocorrelation time for col ", coli, "=", IAcFull
					
				if maxAcor < IAcFull: maxAcor = int(IAcFull)
			#print "Max full autocorr time", maxAcor
		
			# Get moving averages
			self.mvAvgs[seedi] = np.zeros((self.nofDataCols, (self.logData[seedi][0]).size ))
			self.logDataMeanLines = np.ones((self.nofDataCols, (self.logData[seedi][0]).size ))
			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1
				self.logDataMeanLines[dataColi] = self.means[seedi][dataColi]
				if is_odd(maxAcor):
					self.mvAvgs[seedi][dataColi] = moving_average(self.logData[seedi][coli], maxAcor)
				else:
					self.mvAvgs[seedi][dataColi] = moving_average(self.logData[seedi][coli], maxAcor + 1)
		
			# Get equilibration point (moving average vs mean intersection)
			self.eqPoints[seedi] = np.zeros((self.nofDataCols))
			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1
				self.mvAvgs_Means_Xs = intersections(self.mvAvgs[seedi][dataColi], self.logDataMeanLines[dataColi])
				self.mvAvgs_Means_1stX = self.mvAvgs_Means_Xs[0]
				self.eqPoints[seedi][dataColi] = self.mvAvgs_Means_1stX
		
			# Discard equilibration period
			self.eqPoint[seedi] = int(np.max(self.eqPoints[seedi]))
			self.trimLogData[seedi] = self.logData[seedi][:, self.eqPoint[seedi] :]
		#

	def AnalyzeAutocorrelation(self, nofAddMethods, lagMaxima):	
		'''
		Get autocorrelation funtions on production period for
		every data series.
		'''
		for seedi in range(len(self.FNSeeds)):
			self.cuts[seedi] = np.zeros((self.nofDataCols))
			self.Iacs[seedi] = np.zeros((self.nofDataCols))
			self.ESSs[seedi] = np.zeros((self.nofDataCols))
			self.trimAutocorrFunc[seedi] = np.empty((self.nofDataCols, 7, (self.logData[seedi][0]).size)) # new
			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1

				# Autocorrelation function calculated with Grossfield 2009 method
				M = lagMaxima[dataColi] # Calculate up to point M or until it reaches 0
				trimN = self.trimLogData[seedi][coli].size
				if M == 0: M = trimN
	
				autocorrFunc, self.cuts[seedi][dataColi], self.Iacs[seedi][dataColi], self.ESSs[seedi][dataColi] \
					= autocorr1(M, self.trimLogData[seedi][coli])
				self.trimAutocorrFunc[seedi][dataColi][0][:autocorrFunc.size] = autocorrFunc
	
				# Standard error based on autocorrelation time
				#SE = self.stds[seedi][dataColi] * np.sqrt(self.Iacs[seedi][dataColi] / float(trimN))
				#print "Standard error", SE
	
				# Additional methods if required
				add_funcs = [autocorr1, autocorr2, autocorr3, autocorr4, autocorr5]
				for i in range(0, nofAddMethods):
					autocorrFunc = add_funcs[i](M, self.trimLogData[seedi][coli])
					self.trimAutocorrFunc[seedi][dataColi][i][:autocorrFunc.size] = autocorrFunc
	#

	def PyMBARAutocorrelation(self):
		autocor = Autocorrelation.Autocorrelation()
		for seedi in range(len(self.FNSeeds)):

			autocor.getData( (self.logData[seedi][[self.dataCols], :])[0] )
			self.t_g_Neff[seedi] = autocor.pymbarDetectEquilibration_fft()

			dataColi = -1
			for coli in self.dataCols:
				dataColi += 1
			
				IAc = (self.t_g_Neff[seedi][dataColi][1] - 1.0) / 2.0
				print("FFT PyMBAREqPnt statIneff IAc maxESS:")
				print ("eqPointMax(" + str(coli) + ")= " + str(self.t_g_Neff[seedi][dataColi][0]) \
					+ " g= " + str(self.t_g_Neff[seedi][dataColi][1]) \
					+ " IAc= " + str(IAc) \
					+ " Neffmax= " + str(self.t_g_Neff[seedi][dataColi][2]))
	#






