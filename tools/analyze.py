import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from autocorFuncs import *

def is_odd(num):
	return num % 2 != 0
#

def print2D(M):
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			print M[i, j],
		print
#
	
def print2DAndExit(M):
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			print M[i, j],
		print
	exit(0)
#
	
parser = argparse.ArgumentParser()
parser.add_argument('--simDirs', default=None, nargs='+',
	help='Directory with input data files')
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
parser.add_argument('--makeplots', default=[], nargs='+', 
	help='Make plots')
parser.add_argument('--savefig', action='store_true', default=False,
	help='Save the plot into a file')
args = parser.parse_args()

# General plot parameters
colors = ['black', 'red', 'orange', 'blue', 'magenta', 'cyan', 'green', 'pink', 'black', 'black']
if args.makeplots:
	fignum = 0


# Main loop: iterates through seeds = simulation reps
tiny = 0.0000001
for ri in range(len(args.FNSeeds)):
	# Get file names
	FNList = []
	for di in range(len(args.simDirs)):
		FNList.append(glob.glob(os.path.join(args.simDirs[di] + "log." + args.FNSeeds[ri]))[0])
	nfiles = len(FNList)

	# Get raw data
	print "Gathering data from", FNList, "...",
	rawdataChunks = []	
	for li in range(nfiles):
		with open(FNList[li], 'r') as in_FN1:
			rawdataChunks.append(np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, invalid_raise = False))
			rawdataChunks[li] = rawdataChunks[li][::args.stride]
	rawdataChunks = np.array(rawdataChunks)
	rawdata = np.concatenate(rawdataChunks)

	# Transpose data and add another column for acceptance
	logData = np.concatenate((rawdata.transpose(), np.zeros((1, rawdata.shape[0]))), axis=0)
	ncols = logData.shape[0]
	print "Done."

	# Attach acceptance indicator on the last column
	print "Attaching acceptance indicator based on difference of eolumns", 2, 3,
	for i in range(rawdata.shape[0]):
		if rawdata[i][2] != rawdata[i][3]:
			logData[-1][i] = 1
	avgAcc = np.mean(logData[-1])
	print "Done."

	# Save moments for entire chosen timeseries
	nofDataCols = len(args.datacols)
	dataCols = args.datacols

	variances = np.zeros((nofDataCols))
	stds = np.zeros((nofDataCols))
	means = np.zeros((nofDataCols))

	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		print "Column", coli,
		means[dataColi] = np.mean(logData[coli])
		variances[dataColi] = np.var(logData[coli])
		stds[dataColi] = np.std(logData[coli])
		print "mean stdev skew kurt ", means[dataColi], stds[dataColi], scipy.stats.skew(logData[coli]), scipy.stats.kurtosis(logData[coli])

	# Autocorrelation functions for data columns
	maxAcor = 0
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		if args.acf:
			# Autocorrelation time estimate (2009 Grossfield)
			# Compute correlation function for M values or until it reaches 0
			M = args.Ms[dataColi] # Calculate up to point M
			N = logData[coli].size # Sample size
			if M == 0: M = N

			# Autocorrelation functions calculated with different methods
			corrFull, cut = CestGrossfield(M, logData[coli]) # Beyond the cut point it becomes noisy
			
			# Integrated autocorrelation time
			IAcFull = np.sum(corrFull)
			print "Integrated autocorrelation time full", IAcFull
			
			if maxAcor < IAcFull: maxAcor = int(IAcFull)
	print "Max full autocorr time", maxAcor

	# Get moving averages
	maCols = np.zeros((nofDataCols, (logData[0]).size ))
	logDataMeanLines = np.ones((nofDataCols, (logData[0]).size ))
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		logDataMeanLines[dataColi] = means[dataColi]
		if is_odd(maxAcor):
			maCols[dataColi] = moving_average(logData[coli], maxAcor)
		else:
			maCols[dataColi] = moving_average(logData[coli], maxAcor + 1)

	# Get moving average vs mean intersection
	firstIntersect = 0
	diff = prevDiff = 0.0
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		for i in range((logData[0]).size):
			if not np.isnan(maCols[dataColi][i]):
				diff = maCols[dataColi][i] - logDataMeanLines[dataColi][i]
				if (diff * prevDiff) < 0:
					firstIntersect = i
					break
				prevDiff = diff

	# Get equilibration point
	eqPoint = firstIntersect
	print "Equilibration point", eqPoint

	# Get production run (Trim the series)
	trimCols = np.zeros((ncols, (logData[0][eqPoint:]).size ))
	for coli in range(ncols):
		trimCols[coli] = logData[coli][eqPoint:]


	# Get autocorrelation funtions on production runs
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		if args.acf:
			# Autocorrelation time estimate (2009 Grossfield)
			# Compute correlation function for M values or until it reaches 0
			M = args.Ms[dataColi] # Calculate up to point M
			trimN = trimCols[coli].size # Sample size
			if M == 0: M = trimN

			# Autocorrelation functions calculated with different methods
			corr = np.zeros((7, M)) # (# of autocor functions, X size)
			corr[0], cut = CestGrossfield(M, trimCols[coli]) # Beyond the cut point it becomes noisy
			
			# Integrated autocorrelation time
			Iac = np.sum(corr[0])
			ESS = N / Iac
			print "Integrated autocorrelation time", Iac 
			print "Grossfield independent samples", ESS

			# Additional methods
			#add_funcs = [autocorr1, autocorr2, autocorr3, autocorr4, autocorr5]
			add_funcs = [autocorr2]
			for i in range(2, args.nofAddMethods):
				corr[i] = add_funcs[i](M, trimCols[coli])

			# Standard error based on autocorrelation time fitting
			SE = stds[dataColi] * np.sqrt(Iac / float(N))
			print "Standard error", SE

			# Plots
			if(('acf' in args.makeplots) or ('all' in args.makeplots)):
				fignum += 1
				fig = plt.figure(fignum)
				fig.suptitle("Autocorrelation time function") 
				ax = plt.subplot(nofDataCols, 1, dataColi + 1)
				line, = ax.plot(corr[1], label='Col ' +  str(coli) + ' ' + ' f1', color='gray')

				line.set_dashes([2,2,10,2])
				ax.plot(corrFull[0], label='Col ' +  str(coli), color='black')
				ax.plot(corr[0], label='Col ' +  str(coli), color='red')

				ax.legend()
	
				ax.set_xlabel(r'$\mathrm{\tau}$', fontsize=8)
				ax.set_ylabel(r'$\mathrm{C(\tau)}$', fontsize=8)
		  
				plt.setp(ax.get_xticklabels(), rotation='vertical', fontsize=8)
				plt.setp(ax.get_yticklabels(), rotation='horizontal', fontsize=8)
		  
				if args.savefig:
					  figFN = 'temp.acf.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')
	
		if(('data' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			fig = plt.figure(fignum)
			fig.suptitle(FNList[li])
			ax = plt.subplot(nofDataCols, 1, dataColi+1)

			# Initial data
			ax.plot(logData[coli], label='Col ' +  str(coli) + ' real', color='gray')
			# Moving average
			ax.plot(maCols[dataColi], label='Col ' +  str(coli) + ' real', color='green')
			ax.plot(logDataMeanLines[dataColi], label='Col ' +  str(coli) + ' real', color='cyan')
			# Trimmed data
			ax.plot( np.concatenate((np.zeros((eqPoint)), trimCols[coli])), label='Col ' +  str(coli) + ' real', color='black')

			ax.plot(moving_average((1/avgAcc) * means[dataColi] *  logData[-1], maxAcor), color='blue')

			ax.set_xlabel(r'$\mathrm{t}$', fontsize=8)
			ax.set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)

			if args.savefig:
				  figFN = 'temp.data.pdf'
				  plt.savefig(figFN, dpi=600, format='pdf')

if args.makeplots:
	plt.show()












