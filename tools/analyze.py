import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from logAnalyzer import LogAnalyzer
from trajAnalyzer import TrajectoryAnalyzer
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
parser.add_argument('--savefigs', action='store_true', default=False,
	help='Save the plot into a file')
args = parser.parse_args()

# Analyze trajectory
TA = TrajectoryAnalyzer('ala71/ligand.prmtop', 'ala71', args.FNSeeds, args.simDirs)
#TA.ReadPdbs(verbose = False)
TA.ReadDcds(verbose = True)
TA.RMSD()
TA.RG()
TA.SASA()
TA.Helicity()

# General plot parameters
colors = ['red', 'orange', 'magenta', 'cyan', 'green', 'pink']

# Analyze log data	
LA = LogAnalyzer(args.FNSeeds, args.simDirs, args.datacols, args.skip_header, args.skip_footer, args.stride)

# Read files and do basic statistics
LA.Read(verbose = True)

# Find equilibration points
print ("Finding equilibration points...")
LA.FindEquilibrationPoints( lagMaxima = ([0] * len(args.datacols)) )
print ("Done.")

# Compute autocorrelation related quantities
print ("Recalculating autocorrelation functions for production period...")
LA.Autocorrelation(args.nofAddMethods, lagMaxima = ([0] * len(args.datacols)) )
print ("Done.")

# Print
for seedi in range(len(args.FNSeeds)):
	dataColi = -1
	for coli in args.datacols:
		dataColi += 1
		print ("Seed column", args.FNSeeds[seedi], coli),
		#print ("mean stdev skew kurt "),
		#print (LA.means[seedi][dataColi], LA.stds[seedi][dataColi], \
		#	scipy.stats.skew(LA.logData[seedi][coli]), \
		#	scipy.stats.kurtosis(LA.logData[seedi][coli])),

		print ("eqPoint", LA.eqPoints[seedi][dataColi]),
		print ("Iac and ESS", LA.Iacs[seedi][dataColi], LA.ESSs[seedi][dataColi])
	print ("Seed general eqPoint", LA.eqPoint[seedi])

# Make plots
if args.makeplots:
	print ("Plotting...")
	figsPerSeed = 3
	figs = [None] * len(args.FNSeeds) * figsPerSeed
	axes = [None] * len(args.FNSeeds) * figsPerSeed

	# Plot data series
	for seedi in range(len(args.FNSeeds)):
		fignum = -1

		if(('data' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			figIx = (seedi * figsPerSeed) + fignum
			figs[figIx], axes[figIx] = plt.subplots(LA.nofDataCols, 2)
			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				# Plot initial data
				axes[figIx][dataColi, 0].plot(LA.logData[seedi][coli], label='Col ' +  str(coli) + ' real', color='gray')

				# Plot moving averages
				axes[figIx][dataColi, 0].plot(LA.mvAvgs[seedi][dataColi], label='Col ' +  str(coli) + ' real', color='green')
	
				logDataMeanLine = np.ones(((LA.logData[seedi][coli]).size )) * LA.means[seedi][dataColi]
				axes[figIx][dataColi, 0].plot(logDataMeanLine, label='Col ' +  str(coli) + ' real', color='cyan')
	
				# Plot trimmed data
				trim2bPlotted = np.concatenate(( np.full(( LA.eqPoint[seedi]), np.nan)  , LA.trimLogData[seedi][coli] ))
				axes[figIx][dataColi, 0].plot( trim2bPlotted, \
					label='Col ' +  str(coli) + ' real', color='black')

				# Plot acceptance scaled
				scaleFactor = (1/LA.avgAcc[seedi]) * LA.means[seedi][dataColi]
				acc2bPlotted = moving_average(scaleFactor * LA.logData[seedi][-1], int(LA.Iacs[seedi][dataColi]))
				axes[figIx][dataColi, 0].plot(acc2bPlotted, color='blue')
	
				# Format plots
				axes[figIx][dataColi, 0].set_xlabel(r'$\mathrm{t}$', fontsize=8)
				axes[figIx][dataColi, 0].set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)
	
				# Save plots
				if args.savefigs:
					  figFN = 'temp.data.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')
	
		# Plot autocorrelation functions
		if(('acf' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			figIx = (seedi * figsPerSeed) + fignum
			figs[figIx], axes[figIx] = plt.subplots(LA.nofDataCols, 2)
			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				# Plot autocor f on the whole series
				axes[figIx][dataColi, 0].plot(LA.corrFull[seedi][dataColi], \
					label='Col ' +  str(coli) + str(' Grossfield Full'), color='black')

				# Plot autocor f on the production series
				axes[figIx][dataColi, 0].plot(LA.trimCorr[seedi][dataColi][0], \
					label='Col ' +  str(coli) + ' ' + ' Grossfield', color='blue')

				# Plot other autocorrelation functions
				for i in range(0, args.nofAddMethods):
					line, = axes[figIx][dataColi, 0].plot(LA.trimCorr[seedi][dataColi][i], \
						label='Col ' +  str(coli) + ' ' + str(autocorrLabels[i]), color = colors[i])
					line.set_dashes([2,2,10,2])
				
				# Format plots
				axes[figIx][dataColi, 0].legend()
				axes[figIx][dataColi, 0].set_xlabel(r'$\mathrm{\tau}$', fontsize=8)
				axes[figIx][dataColi, 0].set_ylabel(r'$\mathrm{C(\tau)}$', fontsize=8)
				plt.setp(axes[figIx][dataColi, 0].get_xticklabels(), rotation='vertical', fontsize=8)
				plt.setp(axes[figIx][dataColi, 0].get_yticklabels(), rotation='horizontal', fontsize=8)
		  
				# Save plots
				if args.savefigs:
					  figFN = 'temp.acf.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')

		# Plot trajectory based geometric functions
		if(('helicity' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			figIx = (seedi * figsPerSeed) + fignum
			figs[figIx], axes[figIx] = plt.subplots(2, 2)
			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			axes[figIx][0, 0].plot( TA.rmsds[seedi], \
				label='RMSD', color='black')
			axes[figIx][0, 0].legend()

			axes[figIx][0, 1].plot( TA.RGs[seedi], \
				label='RGs', color='black')
			axes[figIx][0, 1].legend()

			axes[figIx][1, 0].plot( TA.helicities1[seedi], \
				label='Helicity', color='black')
			axes[figIx][1, 0].legend()

			axes[figIx][1, 1].scatter( TA.totSASAs[seedi], TA.RGs, \
				label='SASA vs RG', color='black', s = 1**2)
			axes[figIx][1, 1].legend()
if args.makeplots:
	plt.show()
#
#
