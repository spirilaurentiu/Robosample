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

if "traj" in args.analyze:
	# Analyze trajectory
	TA = TrajectoryAnalyzer( (args.molName + '/ligand.prmtop'), args.molName, args.FNSeeds, args.simDirs)
	#TA.ReadPdbs(verbose = False)
	TA.ReadDcds(verbose = True)
	TA.RMSD()
	TA.RG()
	TA.SASA()
	TA.Helicity()

if "log" in args.analyze:
	# Analyze log data	
	LA = LogAnalyzer(args.FNSeeds, args.simDirs, args.datacols, args.skip_header, args.skip_footer, args.stride)
	LA.Read(verbose = True) # Read files and do basic statistics

	# Find equilibration points
	print ("Finding equilibration points...")
	LA.FindEquilibrationPoints( lagMaxima = ([0] * len(args.datacols)) )
	print ("Done.")

	# Compute autocorrelation related quantities
	print ("Recalculating autocorrelation functions for production period...")
	LA.AnalyzeAutocorrelation(args.nofAddMethods, lagMaxima = ([0] * len(args.datacols)) )
	print ("Done.")

	# Compute autocorrelation related quantities
	print ("Recalculating autocorrelation functions with PyMBAR.")
	LA.PyMBARAutocorrelation()
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
	
			print ("eqPointFirst", LA.eqPoints[seedi][dataColi]),
			print ("Iac and ESS", LA.Iacs[seedi][dataColi], LA.ESSs[seedi][dataColi])
		print ("Seed general eqPointFirst", LA.eqPoint[seedi])

fignum = -1
###########################
# SAMPLE STAATISTICS
print('SAMPLE STATISTICS')
###########################

nofSeeds = len(args.FNSeeds)
nofDataCols = len(args.datacols)
## Iterate through data columns
#meanOfMeans = np.empty((nofDataCols))
#stdErrs = np.empty((nofDataCols))
#
#means = np.array(LA.means)
#stds = np.array(LA.stds)
#
## Mean of the sampling distribution of the mean
#meanOfMeans = np.mean(means, axis=0)
#
## Standard error of the mean
#stdErrs = np.std(means, axis=0)
#
#print("Sample means", means)
#print(args.molName, "MeanOfMeans", meanOfMeans[0], meanOfMeans[1])
#print(args.molName, "StandardErrors", stdErrs[0], stdErrs[1])
#
#sampleSize = LA.logData[0][args.datacols[0]].size
#
#dataColi = -1
#for coli in args.datacols:
#	dataColi += 1
#
#	fignum += 1
#	statsFig, statsAxes = plt.subplots(1, 1)
#	plt.suptitle('Stats ' + str(coli))
#
#	meanOfMeansLine = meanOfMeans[dataColi] * np.ones((sampleSize))
#	stdErrsLine = stdErrs[dataColi] * np.ones((sampleSize)) # * np.sqrt(sampleSize)
#	statsAxes.plot(meanOfMeansLine, color='black')
#	line, = statsAxes.plot(meanOfMeansLine + stdErrsLine, color='gray')
#	line.set_dashes([2,2,10,2])
#	line, = statsAxes.plot(meanOfMeansLine - stdErrsLine, color='gray')
#	line.set_dashes([2,2,10,2])
#
#	# Iterate through seeds
#	for seedi in range(nofSeeds):
#		series = LA.logData[seedi][coli]
#		cumAvg = np.empty((series.size))
#		cumStd = np.empty((series.size))
#
#		for i in range(series.size):
#			cumAvg[i] = np.mean(series[0:i+1])
#			cumStd[i] = np.std(series[0:i+1])
#
#		line, = statsAxes.plot(cumAvg)
#		#statsAxes.plot(cumAvg + cumStd, dashes=[6, 2], color = line.get_color())
#		#statsAxes.plot(cumAvg - cumStd, dashes=[6, 2], color = line.get_color())
#		
#plt.show()
#exit(0)

# Make plots
if args.makeplots:
	figsPerSeed = 3
	figs = [None] * len(args.FNSeeds) * figsPerSeed
	axes = [None] * len(args.FNSeeds) * figsPerSeed

	# Plot data series
	for seedi in range(nofSeeds):

		print ('Plotting seed', args.FNSeeds[seedi],  '...')
		eqPoints = {}

		fignum = -1

		if(('log' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			figIx = (seedi * figsPerSeed) + fignum
			figs[figIx], axes[figIx] = plt.subplots(LA.nofDataCols, 3)
			logAxes = axes[figIx]
			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				series = LA.logData[seedi][coli]
				###########################
				# Plot initial data
				###########################

				axes[figIx][dataColi, 0].plot(series, label='Col ' +  str(coli), color='gray')

				# Plot difference quotient
				toBePlotted = difference_quotient(series, args.logDiffWin)
				scaleFactor = np.abs(np.mean(series)) / np.abs(np.mean(toBePlotted)) / 10.0
				toBePlotted = scaleFactor * toBePlotted
				#axes[figIx][dataColi, 0].plot(toBePlotted, label='Col ' +  str(coli) + ' diffQuot', color='pink')

				# Plot diff quatioent X intercept
				XAxisIntersections = intersections(toBePlotted, np.zeros((toBePlotted.size)))
				if XAxisIntersections.size:
					XAxisIntersections = XAxisIntersections - 0
					eqPoints.update({'logcol' + str(coli):  XAxisIntersections[0]})
					print("eqPoint(logcol" + str(coli) + ")", eqPoints['logcol' + str(coli)])
				zeroArray = np.zeros((XAxisIntersections.size))
				axes[figIx][dataColi, 0].scatter(XAxisIntersections, zeroArray, color='red')

				# Plot moving averages
				axes[figIx][dataColi, 0].plot(LA.mvAvgs[seedi][dataColi], label='Col ' +  str(coli) + ' MA', color='pink')

				# Plot the average as a horizontal line	
				logDataMeanLine = np.ones(((LA.logData[seedi][coli]).size )) * LA.means[seedi][dataColi]
				axes[figIx][dataColi, 0].plot(logDataMeanLine, label='Col ' +  str(coli) + ' mean', color='cyan')
	
				# Plot data without the equilibration period
				trim2bPlotted = np.concatenate(( np.full(( LA.eqPoint[seedi]), np.nan)  , LA.trimLogData[seedi][coli] ))
				axes[figIx][dataColi, 0].plot( trim2bPlotted, \
					label='Col ' +  str(coli) + ' trimmed', color='black')

				# Format plots
				axes[figIx][dataColi, 0].set_xlabel(r'$\mathrm{t}$', fontsize=8)
				axes[figIx][dataColi, 0].set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)
				axes[figIx][dataColi, 0].legend()

				###########################
				# Plot acceptance
				###########################
				acceptance = LA.logData[seedi][-1]
				acceptanceMA = moving_average(acceptance, args.accMAWin)


				axes[figIx][dataColi, 1].plot(acceptanceMA, label='Acceptance', color='blue')
				axes[figIx][dataColi, 1].scatter(range(0, acceptance.size), acceptance, s=1, color='black')

				# PLot difference quotient of acceptance
				toBePlotted = difference_quotient(acceptanceMA, args.accDiffWin)
				scaleFactor = np.abs(np.mean(acceptanceMA)) / np.abs(np.mean(toBePlotted))
				scaleFactor = 2000.0
				toBePlotted = scaleFactor * toBePlotted
				#axes[figIx][dataColi, 1].plot(toBePlotted, label='Acc diffQuot', color='pink')

				# Plot diff quatioent X intercept
				XAxisIntersections = intersections(toBePlotted, np.zeros((toBePlotted.size)))
				if XAxisIntersections.size:
					XAxisIntersections = XAxisIntersections - 0
					eqPoints.update({'acc':  XAxisIntersections[0]})
					print("eqPoint(acc)", eqPoints['acc'])
				zeroArray = np.zeros((XAxisIntersections.size))
				axes[figIx][dataColi, 1].scatter(XAxisIntersections, zeroArray, label='Acc MvDiff0s', color='red')
	
				# Format plots
				axes[figIx][dataColi, 1].set_xlabel(r'$\mathrm{t}$', fontsize=8)
				axes[figIx][dataColi, 1].set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)
				axes[figIx][dataColi, 1].legend()
	
				# Save plots
				if args.savefigs:
					  figFN = 'temp.log.pdf'
					  plt.savefig(figFN, dpi=600, format='pdf')
	
				# Detailed analysis of the autocorrelation function
				print('AUTOCORRELATION DETAILED ANALYSIS')
				fignum += 1
				fftFig, fftAxes = plt.subplots(1, 1)
				#X = np.arange(0, 4, 0.1)
				#T = 0.1
				#N = int(4.0 / 0.1)
				#X = np.linspace(0.0, N*T, N)
				#series = np.cos((2 * np.pi) * X) + np.cos(2*X)

				T = 1.0
				N = (LA.logData[seedi][coli][1000:]).size
				X = np.linspace(0.0, N*T, N)
				series = LA.logData[seedi][coli][1000:]
				
				fftAxes.plot(X, series, label='Col ' +  str(coli) + ' ')

				fignum += 1
				fftFig, fftAxes = plt.subplots(1, 1)
				plt.suptitle('FFT ' + str(args.FNSeeds[seedi]))

				(freq, cf, sf, invSf, resultF) = autocorrDetail(0, series)
				#for x in invSf:
				#	print(x.real, x.imag)
				XF = np.linspace(0.0, 1.0/(2.0*T), N/2)
				fftAxes.plot(XF, np.abs(cf[:N//2]), label='Col ' +  str(coli) + ' ')
				fftAxes.set_xlim([0, 0.01])
				#fftAxes.plot(XF, cf.imag[:N//2], label='Col ' +  str(coli) + ' ')

				#x = LA.logData[seedi][coli]
				#result = np.correlate(x, x, mode='full')
				#fftAxes.plot(result)
				#fftAxes.plot(result[result.size // 2:])
				#plt.show()
				##################################################

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
				axes[figIx][dataColi, 0].plot(LA.fullAutocorrFunc[seedi][dataColi], \
					label='Col ' +  str(coli) + str(' fullAutocorrFunc'), color='black')

				# Plot autocor f on the production series
				axes[figIx][dataColi, 0].plot(LA.trimAutocorrFunc[seedi][dataColi][0], \
					label='Col ' +  str(coli) + ' ' + ' trimAutocorrFunc', color='blue')

				# Plot other autocorrelation functions
				for i in range(0, args.nofAddMethods):
					line, = axes[figIx][dataColi, 0].plot(LA.trimAutocorrFunc[seedi][dataColi][i], \
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
		if(('traj' in args.makeplots) or ('all' in args.makeplots)):
			fignum += 1
			figIx = (seedi * figsPerSeed) + fignum
			figs[figIx], axes[figIx] = plt.subplots(2, 2)
			plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			series = [TA.rmsds[seedi], TA.RGs[seedi], TA.helicities1[seedi], TA.totSASAs[seedi]]
			seriesLabels = ['RMSD', 'RG', 'Helicity', 'SASA']
			for si in range(len(series)):
				axes[figIx][int(si/2), (si%2)].plot(series[si], label=seriesLabels[si], color='black')

				# PLot difference quotient of acceptance
				toBePlotted = difference_quotient(series[si], args.trajDiffWin)
				scaleFactor = np.abs(np.mean(series[si])) / np.abs(np.mean(toBePlotted)) / 30.0
				toBePlotted = scaleFactor * toBePlotted
				#axes[figIx][int(si/2), (si%2)].plot(toBePlotted, label=seriesLabels[si] + 'diffQuot', color='pink')

				# Plot diff quatioent X intercept
				XAxisIntersections = intersections(toBePlotted, np.zeros((toBePlotted.size)))
				if XAxisIntersections.size:
					XAxisIntersections = XAxisIntersections - 0
					print("eqPoint(" + seriesLabels[si] + ")", XAxisIntersections[0])
				zeroArray = np.zeros((XAxisIntersections.size))
				axes[figIx][int(si/2), (si%2)].scatter(XAxisIntersections, zeroArray, label=seriesLabels[si] + 'diffQuot0s', color='red')
				axes[figIx][int(si/2), (si%2)].legend()

			#axes[figIx][1, 1].scatter( TA.totSASAs[seedi], TA.RGs, \
			#	label='SASA vs RG', s = 1**2, cmap='PuBu_r')
			#axes[figIx][1, 1].legend()

		PyMBARXAxisIntersections = len(args.datacols) * [None]
		if(('pymbar' in args.makeplots) or ('all' in args.makeplots)):
			#fignum += 1
			#figIx = (seedi * figsPerSeed) + fignum
			#figs[figIx], axes[figIx] = plt.subplots(2, 2)
			#plt.suptitle('Seed ' + str(args.FNSeeds[seedi]))

			dataColi = -1
			for coli in args.datacols:
				dataColi += 1

				series = LA.t_g_Neff[seedi][dataColi][3]
				#axes[figIx][dataColi, 0].set_ylim((np.min(series) - np.std(series), np.max(series) + np.std(series)))
				#axes[figIx][dataColi, 0].plot(series, label='Col ' +  str(coli), color='black')
				logAxes[dataColi, 2].set_ylim((np.min(series) - np.std(series), np.max(series) + np.std(series)))
				logAxes[dataColi, 2].plot(series, label='Col ' +  str(coli), color='black')
	
				# Plot difference quotient
				toBePlotted = difference_quotient(series, 301)
				scaleFactor = np.abs(np.mean(series)) / np.abs(np.mean(toBePlotted)) / 100.0
				toBePlotted = scaleFactor * toBePlotted
				#logAxes[dataColi, 2].plot(toBePlotted, label='Col ' +  str(coli) + ' diffQuot', color='pink')

				# Plot diff quatioent X intercept
				PyMBARXAxisIntersections[dataColi] = intersections(toBePlotted, np.zeros((toBePlotted.size)))
				if PyMBARXAxisIntersections[dataColi].size:
					PyMBARXAxisIntersections[dataColi] = PyMBARXAxisIntersections[dataColi] - 0
					eqPoints.update({'pymbar' + str(coli):  PyMBARXAxisIntersections[dataColi][0]})
					print("eqPoint(pymbar" + str(coli) + ")", eqPoints['pymbar' + str(coli)])
				zeroArray = np.zeros((PyMBARXAxisIntersections[dataColi].size))
				logAxes[dataColi, 2].scatter(PyMBARXAxisIntersections[dataColi], zeroArray, color='red')
				logAxes[dataColi, 2].legend()


				# Calc 
				colIndex = dataColi
				print("Calculating for column with index", colIndex)
				key_max = max(eqPoints.keys(), key=(lambda k: eqPoints[k]))
				#finalEqPoint = eqPoints[key_max]
				finalEqPoint = eqPoints['logcol' + str(2)]
		
				# Search for certain configurations
				confSeriesEE = LA.logData[seedi][args.datacols[1]] # E-E dist
				confEqPointEE = confSeriesEE.size - 1
				for i in range(confSeriesEE.size):
					if confSeriesEE[i] < 2.0:
						confEqPointEE = i
						break
				print("EE dist went below 2.0 at", confEqPointEE)
		
				if confEqPointEE == confSeriesEE.size - 1:
					print("Simulation not equilibrated.")
		
				# Search for maxima in Neff function around point confEqPoint method 1
				ESSfunction = LA.t_g_Neff[seedi][colIndex][3]
				peaks, _ = scisig.find_peaks(ESSfunction, distance=1000)

				# TODO: move after eq point and work on it
				valleys, _ = scisig.find_peaks(-1.0 * ESSfunction, distance=1000)
		
				finalEqPoint = ESSfunction.size - 1
				# Find stable points to the left then to the right
				for i in range(confEqPointEE, 0, -1): # Search to the left
					if i in PyMBARXAxisIntersections[0]:
						#print("ESSfunction leveled at", finalEqPoint, "to the left")
						pass
					if i in peaks:
						finalEqPoint = i
						print("ESSfunction peaked at", finalEqPoint, "to the left")
						break
		
				# No level to the left
				if finalEqPoint == (ESSfunction.size - 1):
					for i in range(confEqPointEE, ESSfunction.size): # Search to the right
						if i in PyMBARXAxisIntersections[0]:
							#print("ESSfunction leveled at", finalEqPoint, "to the right")
							pass
						if i in peaks:
							finalEqPoint = i
							print("ESSfunction peaked at", finalEqPoint, "to the right")
							break

				# Method 3
				finalEqPoint = np.argmax(ESSfunction[0:-1000])
	
				print("Final equilibration point seed",  args.FNSeeds[seedi], "col", coli\
					, finalEqPoint)
				print("Final autocorrelation time seed", args.FNSeeds[seedi], "col", coli\
					, (ESSfunction.size - finalEqPoint) / ESSfunction[finalEqPoint])
	
				trimESSfunction = np.empty((ESSfunction.size)) * np.nan # 0 is the potential energy index
				trimESSfunction[finalEqPoint:] = ESSfunction[finalEqPoint:]
	
				#print1DNPArray(ESSfunction[:300])
				#print1DNPArray(trimESSfunction[:300])
	
				logAxes[dataColi, 2].plot(range(0, ESSfunction.size), ESSfunction, color='black')
				logAxes[dataColi, 2].plot(range(0, ESSfunction.size), trimESSfunction, color='orange')
				logAxes[dataColi, 2].plot(peaks, ESSfunction[peaks], "x")
				logAxes[dataColi, 2].plot(valleys, ESSfunction[valleys], "x")

if args.makeplots:
	plt.show()
#
#
