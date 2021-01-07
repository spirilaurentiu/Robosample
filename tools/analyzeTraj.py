# Robosample tools
import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse

#from scipy.signal import find_peaks
import scipy.signal as scisig

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from logAnalyzer import LogAnalyzer
from trajAnalyzer import TrajectoryAnalyzer
from autocorFuncs import *
from jumpDetect import *

def print1DNPArray(series):
	print("\n")
	for k in range(series.shape[0]):
		print(series[k], end = ' ')
		if k % 80 == 79:
			print("\n")
	print("\n")
#

parser = argparse.ArgumentParser()
parser.add_argument('--simDirs', default=None, nargs='+',
	help='Directory with input data files')
parser.add_argument('--logDirs', default=None, nargs='+',
	help='Directory with input data files')
parser.add_argument('--molName', default=None, 
	help='Molecule name also the name of the directory of the prmtop file.')
parser.add_argument('--FNSeeds', default=None, nargs='+',
	help='Name root for first input data files')
parser.add_argument('--skip_header', default=0, type=int,
	help='# of lines to be skipped by numpy.loadtxt from the beggining')
parser.add_argument('--skip_footer', default=0, type=int,
	help='# of lines to be skipped by numpy.loadtxt at the end')
parser.add_argument('--datacols', default=0, type=int, nargs='+',
	help='Columns to be read by numpy.loadtxt')
parser.add_argument('--stride', default=1, type=int,
	help='Stride for the read lines.')
parser.add_argument('--nbins', default=50, type=int,
	help='Number of bins for histograms. Default for 7.2 spacing')
parser.add_argument('--analyze', default=[], nargs='+', 
	help='Decide what to analyze')
parser.add_argument('--Ms', default=[0], type=int, nargs='+', 
	help='Number of steps to sum the correlation function on. One M per column')
parser.add_argument('--fitacf', action='store_true', default=False,
	help='Fit autocorrelation function to an exponential')
parser.add_argument('--acf', action='store_true', default=False,
	help='Compute autocorrelation function.')
parser.add_argument('--bse', action='store_true', default=False,
	help='Use block averaging')
parser.add_argument('--nofAddMethods', default=0, type=int,
	help='Number of additional methods wrote by James on StackOverflow.')
parser.add_argument('--logDiffWin', default=3, type=int,
	help='Moving difference window.')
parser.add_argument('--trajDiffWin', default=3, type=int,
	help='Moving difference window.')
parser.add_argument('--accMAWin', default=3, type=int,
	help='Moving average window.')
parser.add_argument('--accDiffWin', default=3, type=int,
	help='Moving difference window.')
parser.add_argument('--makeplots', default=[], nargs='+', 
	help='Make plots')
parser.add_argument('--savefigs', action='store_true', default=False,
	help='Save the plot into a file')
args = parser.parse_args()
# General plot parameters
colors = ['red', 'orange', 'magenta', 'cyan', 'green', 'pink']

#from pymbar import timeseries as ts
#from pymbar.utils import ParameterError

nofSeeds = len(args.FNSeeds)
if nofSeeds % 3:
	print(args.FNSeeds)
	print("Nof seeds must be divisible by 3. Exiting...")
	exit(1)

nofParamSets = nofSeeds // 3

seeds = np.array(np.array_split(np.array([int(seed) for seed in args.FNSeeds]) , nofParamSets))
print(seeds)
exit(0)

mofm = np.empty((3990))
sofm = np.empty((3990))
TA = []
fignum = -1
if "traj" in args.analyze:
	for diri in range(len(args.simDirs)):
		# Analyze trajectory
		TA.append(TrajectoryAnalyzer( (args.molName + '/ligand.prmtop'), args.molName, args.FNSeeds, [args.simDirs[diri]] ))
		#TA.ReadPdbs(verbose = False)
		TA[diri].ReadDcds(verbose = True)
		TA[diri].Distance(np.array([[10, 6]]))
	
	#	TA.RMSD()
	#	TA.RG()
	#	TA.SASA()
	#	TA.Helicity()
	for paramSeti in nofParamSets:	
		for framei in range(0, 3990):
			distances = []
			for seedi in range(seeds[paramSeti].shape[0]):
				for diri in range(len(args.simDirs)):
					#print(args.simDirs[diri])
					#print(seedi)
					#print(TA[diri].distances[seedi].flatten()[0:framei])
					distances.append(TA[diri].distances[seedi].flatten()[0:framei+1])
			means = [np.mean(sim) for sim in distances]
			mofm[framei] = np.mean(means)
			sofm[framei] = np.std(means)
		
		plt.plot(mofm)

		
		# Plot trajectory based geometric functions
	#	figsPerSeed = 3
	#	figs = [None] * len(args.FNSeeds) * figsPerSeed
	#	axes = [None] * len(args.FNSeeds) * figsPerSeed
	#	for seedi in range(nofSeeds):
	#		if(('traj' in args.makeplots) or ('all' in args.makeplots)):
	#			fignum += 1
	#			figIx = (seedi * 2) + fignum
	#			figs[figIx], axes[figIx] = plt.subplots(2, 1)
	#			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))
	#	
	#			series = [TA.distances[seedi]]
	#			seriesLabels = ['Distance']
	#			for si in range(len(series)):
	#				axes[figIx][si].plot(series[si], label=seriesLabels[si], color='black')
				
	#			# 
	#			fignum += 1
	#			figIx = (seedi * figsPerSeed) + fignum
	#			figs[figIx], axes[figIx] = plt.subplots(2, 2)
	#			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))
	#	
	#			series = [TA.rmsds[seedi], TA.RGs[seedi], TA.helicities1[seedi], TA.totSASAs[seedi]]
	#			seriesLabels = ['RMSD', 'RG', 'Helicity', 'SASA']
	#			for si in range(len(series)):
	#				axes[figIx][int(si/2), (si%2)].plot(series[si], label=seriesLabels[si], color='black')
	#	
	#				# PLot difference quotient of acceptance
	#				toBePlotted = difference_quotient(series[si], args.trajDiffWin)
	#				scaleFactor = np.abs(np.mean(series[si])) / np.abs(np.mean(toBePlotted)) / 30.0
	#				toBePlotted = scaleFactor * toBePlotted
	#				#axes[figIx][int(si/2), (si%2)].plot(toBePlotted, label=seriesLabels[si] + 'diffQuot', color='pink')
	#	
	#				# Plot diff quatioent X intercept
	#				XAxisIntersections = intersections(toBePlotted, np.zeros((toBePlotted.size)))
	#				if XAxisIntersections.size:
	#					XAxisIntersections = XAxisIntersections - 0
	#					print("eqPoint(" + seriesLabels[si] + ")", XAxisIntersections[0])
	#				zeroArray = np.zeros((XAxisIntersections.size))
	#				axes[figIx][int(si/2), (si%2)].scatter(XAxisIntersections, zeroArray, label=seriesLabels[si] + 'diffQuot0s', color='red')
	#				axes[figIx][int(si/2), (si%2)].legend()
	#	
	#			#axes[figIx][1, 1].scatter( TA.totSASAs[seedi], TA.RGs, \
	#			#	label='SASA vs RG', s = 1**2, cmap='PuBu_r')
	#			#axes[figIx][1, 1].legend()



if args.makeplots:
	plt.show()
#
#
