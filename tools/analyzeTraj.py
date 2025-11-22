# Robosample tools
import sys, os, glob
import numpy as np
import copy
import scipy 
import scipy.stats
import argparse

#from scipy.signal import find_peaks
import scipy.signal as scisig
from scipy.optimize import curve_fit

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from logAnalyzer import LogAnalyzer
from trajAnalyzer import TrajectoryAnalyzer
from autocorFuncs import *
from jumpDetect import *

# Function to print big arrays the way I want
def print1DNPArray(series):
	print("\n")
	for k in range(series.shape[0]):
		print(series[k], end = ' ')
		if k % 80 == 79:
			print("\n")
	print("\n")
#

# Fitting function
def func(x, a, b, c):
    return a * np.exp(-b * x) + c
#

# Arguments
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

parser.add_argument('--distance', default=[0, 1], type=int, nargs='+',
	help='Atom indeces to compute distance for.')

parser.add_argument('--trajDiffWin', default=3, type=int,
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

# [[s1 s2 s3], [s4 s5 6]]
seeds = np.array(np.array_split(np.array([int(seed) for seed in args.FNSeeds]) , nofParamSets))
print(seeds)

cumulativeMean = np.empty((3990))
cumulativeVar = np.empty((3990))
TA = []
fignum = -1
if "traj" in args.analyze:
	for diri in range(len(args.simDirs)):
		# Analyze trajectory
		TA.append(TrajectoryAnalyzer( (args.molName + '/ligand.prmtop'), args.molName, args.FNSeeds, [args.simDirs[diri]] ))
		#TA.ReadPdbs(verbose = False)
		TA[diri].ReadDcds(verbose = True)
		TA[diri].Distance(np.array([[args.distance[0], args.distance[1]]]))
	
	#	TA.RMSD() TA.RG() TA.SASA() TA.Helicity()

	# Put data in Numpy array: dim1 is dir, dim2 is seed, dim3 is data
	# seed 0 dir 0 --------
	# seed 0 dir 1 --------
	# seed 1 dir 0 --------
	# seed 1 dir 1 --------
	# seed 2 dir 0 --------
	# seed 2 dir 1 --------
	data = []
	dim1ix = -1
	for seedi in range(len(args.FNSeeds)):
		for diri in range(len(args.simDirs)):
			dim1ix += 1
			#print("dir", diri, args.simDirs[diri])
			#print("seedi to seed", seedi, args.FNSeeds[seedi])
			print("append TA[", args.simDirs[diri], "].data[", args.FNSeeds[seedi], "]")
			#print(TA[diri].data[seedi].flatten())
			data.append(TA[diri].distances[seedi].flatten())

	# Get maximum length
	lens = [len(datum) for datum in data]
	
	# Put data in a Numpy array with maximum dimensions
	npdata = np.empty((len(data), max(lens))) * np.nan
	for i in range(len(data)):
		for j in range(len(data[i])):
			npdata[i, j] = data[i][j]
	data = npdata
	print("data shape", data.shape)
	print("data")
	print(data)
	totalNofSims = int(data.shape[0])

	# Calculate
	pltNofRows = 2
	pltNofCols = 1

	nofSimsPerParamSet = int(data.shape[0] / nofParamSets)
	print("There are", nofSimsPerParamSet, "simulations per parameter set.")
	cumulativeBias = np.empty(data.shape)
	cumulativeMean = np.empty(data.shape)
	cumulativeVar = np.empty(data.shape)
	print("nofParamSets", nofParamSets)
	print("data.shape", data.shape)
	bias = np.empty((nofParamSets, data.shape[1]))
	sofb = np.empty((nofParamSets, data.shape[1]))
	mofm = np.empty((nofParamSets, data.shape[1]))
	vofm = np.empty((nofParamSets, data.shape[1]))
	sofm = np.empty((nofParamSets, data.shape[1]))
	MSE = np.empty((nofParamSets, data.shape[1]))

	# Get the true value estimate (last average of the longest simulations)
	# Get the latest of each type
	lastValues = []
	for datumi in range(data.shape[0]):
		if not np.isnan(data[datumi][-1]):
			lastValues.append(data[datumi][-1])
	trueValueEst = np.mean(lastValues)
	print("trueValueEst", trueValueEst)

	# Get cumulative average bias and variance
	fig, axs = plt.subplots(pltNofRows, pltNofCols)
	for paramSeti in range(nofParamSets):
		startSim = paramSeti * nofSimsPerParamSet
		print("Parameter set ", paramSeti, "with seeds", seeds[paramSeti])
		for framei in range(0, data.shape[1]):
			paramSetSimsRange = range(startSim, nofSimsPerParamSet + startSim)
			for simi in paramSetSimsRange:
				# Make sure the current frame is not empty
				if data[simi, framei] != None:
					cumulativeMean[simi, framei] = np.mean(data[simi, 0:framei+1])
					cumulativeBias[simi, framei] = (cumulativeMean[simi, framei] - trueValueEst)
					cumulativeVar[ simi, framei] = np.var(data[simi, 0:framei+1])
		
			#print("cumulativeMean.shape", cumulativeMean.shape)
			# Choose includedSims
			#includedSims = paramSetSimsRange
			includedSims = []
			for simi in paramSetSimsRange:
				if data[simi, framei] != None:
					includedSims.append(simi)
			mofm[paramSeti, framei] = np.mean(cumulativeMean[includedSims, framei])
			vofm[paramSeti, framei] = np.var(cumulativeMean[includedSims, framei])
			bias[paramSeti, framei] = mofm[paramSeti, framei] - trueValueEst
			sofm[paramSeti, framei] = np.std(cumulativeMean[includedSims, framei])
			MSE[paramSeti, framei]  = vofm[paramSeti, framei] + (bias[paramSeti, framei])**2

			sofb[paramSeti, framei] = np.std(cumulativeBias[:, framei])
			

	# Compute the half life point of MSE
	halfLifeMSE = np.zeros((nofParamSets), dtype=int)
	for paramSeti in range(nofParamSets):
		hly = np.mean([np.nanmin(MSE[paramSeti]), np.nanmax(MSE[paramSeti])])
		for framei in range(1, data.shape[1] - 1):
			prevDelta = MSE[paramSeti, framei-1] - hly
			nextDelta = MSE[paramSeti, framei+1] - hly 
			if np.sign(prevDelta * nextDelta) == -1:
				halfLifeMSE[paramSeti] = framei
				break

	# Exponential fitting
	# First get the nans out
	noNanMSE = MSE # shallow copy
	for paramSeti in range(nofParamSets):
		for framei in range(noNanMSE.shape[1] - 1):
			if np.isnan(noNanMSE[paramSeti, framei + 1]):
				noNanMSE[paramSeti, framei+1:] = noNanMSE[paramSeti, framei]
				break
	
	# Fit a polynomial of degree 1 to the log
	expFitFactors = np.zeros((nofParamSets, 2), dtype=float)
	for paramSeti in range(nofParamSets):
		expFitFactors[paramSeti] = np.polyfit( np.array(range(noNanMSE.shape[1])), np.log(noNanMSE[paramSeti]), 1)
		#print(expFitFactors[paramSeti])

	# Use fitting functions
	poptpcov = [None] * nofParamSets
	for paramSeti in range(nofParamSets):
		popt, pcov = curve_fit(func, np.array(range(1, noNanMSE.shape[1])), noNanMSE[paramSeti][1:])
		poptpcov[paramSeti] = []
		poptpcov[paramSeti].append(popt)
		poptpcov[paramSeti].append(pcov)
	print(args.molName, args.FNSeeds, end = ' ')
	print("halfLifeMSE expFactors ratios second/first", halfLifeMSE[1]/halfLifeMSE[0], (-1.0 * poptpcov[1][0][1]) / (-1.0 * poptpcov[0][0][1]))

	# Plots
	if args.makeplots:
		print("Plotting cumulative means and stds")
		colors = ['black', 'red']
		for paramSeti in range(nofParamSets):
			startSim = paramSeti * nofSimsPerParamSet
			X = range(0, cumulativeMean.shape[1])
			for simi in range(startSim, nofSimsPerParamSet + startSim):
	
				#X = range(0, cumulativeMean.shape[1])
				Y = cumulativeMean[simi]
				axs[0].plot(X, Y, color = colors[paramSeti])
				axs[0].set_title('Running Means ' + args.molName + ' seed ' + str(args.FNSeeds[0]))
				#axs[0].set_xlim(0, 4000)
				#axs[0].set_ylim(0, 24)
	
				#Y = data[simi]
	
			Y = noNanMSE[paramSeti]
			Yerr = sofb[paramSeti] / np.sqrt(nofSimsPerParamSet)
			#axs[1].plot(X, Y, color = colors[paramSeti])
			axs[1].errorbar(X, Y, yerr=Yerr, errorevery=100, color = colors[paramSeti])
			axs[1].set_title('MSE')
			#axs[1].set_xlim(0, 4000)
			#axs[1].set_ylim(0, 400)

			# Plot the exponential decay
			a = expFitFactors[paramSeti][0]
			b = expFitFactors[paramSeti][1]
			Y = np.exp(a * X) * np.exp(b)
			#axs[1].plot(X, Y, color = colors[paramSeti], linestyle='dashed')

			# Plot the fitted function
			axs[1].plot(X, func(X, *(poptpcov[paramSeti][0])), 'r-', label="Fitted Curve", color = colors[paramSeti], linestyle='dashed')
	
if args.makeplots:
	plt.legend()
	plt.show()
#
#
