import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from autocorFuncs import *

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
colors = ['red', 'orange', 'magenta', 'cyan', 'green', 'pink']
if args.makeplots:
	fignum = 0


# Main loop: iterates through seeds = simulation reps
tiny = 0.0000001

#figs = [None] * len(args.FNSeeds)
#axes = [None] * len(args.FNSeeds)

logData = [None] * len(args.FNSeeds)
trimLogData = [None] * len(args.FNSeeds)
means = [None] * len(args.FNSeeds)
variances = [None] * len(args.FNSeeds)
stds = [None] * len(args.FNSeeds)
mvAvgs = [None] * len(args.FNSeeds)
avgAcc = [None] * len(args.FNSeeds)
eqPoints = [None] * len(args.FNSeeds)
corrFull = [None] * len(args.FNSeeds)
trimCorr = [None] * len(args.FNSeeds)

nofDataCols = len(args.datacols)
dataCols = args.datacols

for seedi in range(len(args.FNSeeds)):
	# Get file names
	FNList = []
	for di in range(len(args.simDirs)):
		FNList.append(glob.glob(os.path.join(args.simDirs[di] + "log." + args.FNSeeds[seedi]))[0])
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
	logData[seedi] = np.concatenate((rawdata.transpose(), np.zeros((1, rawdata.shape[0]))), axis=0)
	ncols = logData[seedi].shape[0]
	print "Done."

	# Attach acceptance indicator on the last column
	print "Attaching acceptance indicator based on difference of eolumns", 2, 3,
	for i in range(rawdata.shape[0]):
		if rawdata[i][2] != rawdata[i][3]:
			logData[seedi][-1][i] = 1
	avgAcc[seedi] = np.mean(logData[seedi][-1])
	print "Done."

	# Save moments for entire chosen timeseries
	variances[seedi] = np.zeros((nofDataCols))
	stds[seedi] = np.zeros((nofDataCols))
	means[seedi] = np.zeros((nofDataCols))
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		print "Column", coli,
		means[seedi][dataColi] = np.mean(logData[seedi][coli])
		variances[seedi][dataColi] = np.var(logData[seedi][coli])
		stds[seedi][dataColi] = np.std(logData[seedi][coli])
		print "mean stdev skew kurt ", means[seedi][dataColi], stds[seedi][dataColi], \
			scipy.stats.skew(logData[seedi][coli]), scipy.stats.kurtosis(logData[seedi][coli])

	# Autocorrelation functions for data columns
	print "Calculating autocorrelation functions..."
	corrFull[seedi] = np.empty((nofDataCols, (logData[seedi][0]).size))
	maxAcor = 0
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		if args.acf:
			# Autocorrelation time estimate (2009 Grossfield)
			M = args.Ms[dataColi] # Calculate up to point M or until it reaches 0
			N = logData[seedi][coli].size # Sample size
			if M == 0: M = N

			# Autocorrelation functions calculated with different methods
			corrFull[seedi][dataColi], cut, IAcFull, ESSfull = CestGrossfield(M, logData[seedi][coli]) 
			print "Full autocorrelation time for col ", coli, "=", IAcFull
			
			if maxAcor < IAcFull: maxAcor = int(IAcFull)
	print "Done."
	print "Max full autocorr time", maxAcor

	# Get moving averages
	mvAvgs[seedi] = np.zeros((nofDataCols, (logData[seedi][0]).size ))
	logDataMeanLines = np.ones((nofDataCols, (logData[seedi][0]).size ))
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		logDataMeanLines[dataColi] = means[seedi][dataColi]
		if is_odd(maxAcor):
			mvAvgs[seedi][dataColi] = moving_average(logData[seedi][coli], maxAcor)
		else:
			mvAvgs[seedi][dataColi] = moving_average(logData[seedi][coli], maxAcor + 1)

	# Get equilibration point (moving average vs mean intersection)
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		mvAvgs_Means_Xs = intersections(mvAvgs[seedi][dataColi], logDataMeanLines[dataColi])
	mvAvgs_Means_1stX = mvAvgs_Means_Xs[0]
	eqPoints[seedi] = mvAvgs_Means_1stX
	print "Equilibration point", eqPoints[seedi]

	# Discard equilibration period
	trimLogData[seedi] = logData[seedi][:, eqPoints[seedi]:]

	# Get autocorrelation funtions on production period
	print "Recalculating autocorrelation functions for production period..."
	trimCorr[seedi] = np.empty((nofDataCols, 7, (logData[seedi][0]).size)) # new
	dataColi = -1
	for coli in dataCols:
		dataColi += 1
		if args.acf:
			# Autocorrelation time estimate (2009 Grossfield)
			M = args.Ms[dataColi] # Calculate up to point M or until it reaches 0
			trimN = trimLogData[seedi][coli].size
			if M == 0: M = trimN

			# Autocorrelation function calculated with Crossfield method
			#trimCorr[seedi][dataColi] = np.zeros((7, trimN)) # (# of autocor functions, X size)
			x, cut, Iac, ESS = CestGrossfield(M, trimLogData[seedi][coli]) # Beyond the cut point it becomes noisy
			trimCorr[seedi][dataColi][0][:x.size] = x
			print "Autocorrelation time and ESS", Iac, ESS

			# Standard error based on autocorrelation time
			SE = stds[seedi][dataColi] * np.sqrt(Iac / float(trimN))
			print "Standard error", SE

			# Additional methods if required
			add_funcs = [autocorr1, autocorr2, autocorr3, autocorr4, autocorr5]
			for i in range(0, args.nofAddMethods):
				x = add_funcs[i](M, trimLogData[seedi][coli])
				trimCorr[seedi][dataColi][i][:x.size] = x


# Make plots
if args.makeplots:
	figs = [None] * len(args.FNSeeds)
	axes = [None] * len(args.FNSeeds)

	# Plot data series
	for seedi in range(len(args.FNSeeds)):
		figs[seedi], axes[seedi] = plt.subplots(nofDataCols, 2)
		plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

		if(('data' in args.makeplots) or ('all' in args.makeplots)):
			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				# Plot initial data
				axes[seedi][dataColi, 0].plot(logData[seedi][coli], label='Col ' +  str(coli) + ' real', color='gray')
				# Plot moving averages
				axes[seedi][dataColi, 0].plot(mvAvgs[seedi][dataColi], label='Col ' +  str(coli) + ' real', color='green')
	
				logDataMeanLine = np.ones(((logData[seedi][coli]).size )) * means[seedi][dataColi]
				axes[seedi][dataColi, 0].plot(logDataMeanLine, label='Col ' +  str(coli) + ' real', color='cyan')
	
				# Plot trimmed data
				axes[seedi][dataColi, 0].plot( np.concatenate((np.zeros((eqPoints[seedi])), trimLogData[seedi][coli])), label='Col ' +  str(coli) + ' real', color='black')
				# Plot acceptance scaled
				scaleFactor = (1/avgAcc[seedi]) * means[seedi][dataColi]
				acc2bPlotted = moving_average(scaleFactor * logData[seedi][-1], maxAcor)
				axes[seedi][dataColi, 0].plot(acc2bPlotted, color='blue')
	
				# Format plots
				axes[seedi][dataColi, 0].set_xlabel(r'$\mathrm{t}$', fontsize=8)
				axes[seedi][dataColi, 0].set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)
	
				# Save plots
				if args.savefig:
					  figFN = 'temp.data.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')
	
		# Plot autocorrelation functions
		if(('acf' in args.makeplots) or ('all' in args.makeplots)):
			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				# Plot autocor f on the whole series
				axes[seedi][dataColi, 0].plot(corrFull[seedi][dataColi], label='Col ' +  str(coli) + str(' Grossfield Full'), color='black')
				# Plot autocor f on the production series
				axes[seedi][dataColi, 0].plot(trimCorr[seedi][dataColi][0], label='Col ' +  str(coli) + ' ' + ' Grossfield', color='blue')
				# Plot other autocorrelation functions
				for i in range(0, args.nofAddMethods):
					line, = axes[seedi][dataColi, 0].plot(trimCorr[seedi][dataColi][i], \
						label='Col ' +  str(coli) + ' ' + str(autocorrLabels[i]), color = colors[i])
					line.set_dashes([2,2,10,2])
				
				# Format plots
				axes[seedi][dataColi, 0].legend()
				axes[seedi][dataColi, 0].set_xlabel(r'$\mathrm{\tau}$', fontsize=8)
				axes[seedi][dataColi, 0].set_ylabel(r'$\mathrm{C(\tau)}$', fontsize=8)
				plt.setp(axes[seedi][dataColi, 0].get_xticklabels(), rotation='vertical', fontsize=8)
				plt.setp(axes[seedi][dataColi, 0].get_yticklabels(), rotation='horizontal', fontsize=8)
		  
				# Save plots
				if args.savefig:
					  figFN = 'temp.acf.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')
if args.makeplots:
	plt.show()












