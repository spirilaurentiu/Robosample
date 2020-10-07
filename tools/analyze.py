import sys, os, glob
import numpy as np
import scipy 
import scipy.stats
import argparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from autocorFuncs import *

parser = argparse.ArgumentParser()
parser.add_argument('--dir', default=None,
	help='Directory with input data files')
parser.add_argument('--inFNRoots', default=None, nargs='+',
	help='Name root for first input data files')
parser.add_argument('--skip_header', default=0, type=int,
	help='# of lines to be skipped by numpy.loadtxt from the beggining')
parser.add_argument('--skip_footer', default=0, type=int,
	help='# of lines to be skipped by numpy.loadtxt at the end')
parser.add_argument('--usecols', default=0, type=int, nargs='+',
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

def is_odd(num):
	return num % 2 != 0
#

# General plot parameters
colors = ['black', 'red', 'orange', 'blue', 'magenta', 'cyan', 'green', 'pink', 'black', 'black']
if args.makeplots:
	fignum = 0


# Main loop
tiny = 0.0000001
for ri in range(len(args.inFNRoots)): # Iterate through roots
	FNlist = glob.glob(os.path.join(args.dir, args.inFNRoots[ri]))
	nfiles = len(FNlist)
	hists = np.zeros((nfiles, 2, args.nbins))
	relHists = np.zeros((nfiles, 2, args.nbins))
	#print FNlist
	
	for li in range(nfiles):
		with open(FNlist[li], 'r') as in_FN1:
			print "file=", FNlist[li]
			#alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, usecols=args.usecols, invalid_raise = False )
			alldata = np.genfromtxt(in_FN1, skip_header=args.skip_header, skip_footer=args.skip_footer, invalid_raise = False )
			alldata = alldata[::args.stride]

			# Reshape data (transpose)
			if np.size(alldata.shape) == 1:
				cols = np.zeros((1, alldata.shape[0]))
				cols[0] = alldata
			else:
				shape_list = list(alldata.transpose().shape)
				shape_list[0] = shape_list[0] + 1 # 
				shape_tuple = tuple(shape_list)
				cols = np.zeros(shape_tuple)
				for coli in range((alldata.transpose()).shape[0]):
					cols[coli] = alldata.transpose()[coli]
			ncols = cols.shape[0]

			# Attach acceptance indicator on the last column
			for i in range(alldata.shape[0]):
				if alldata[i][2] != alldata[i][3]:
					cols[-1][i] = 1
			avgAcc = np.mean(cols[-1])

			#for i in range(cols.shape[0]):
			#	for j in range(cols.shape[1]):
			#		print cols[i, j],
			#	print

			# Save moments for entire timeseries
			nofDataCols = len(args.usecols)
			dataCols = args.usecols

			variances = np.zeros((nofDataCols))
			stds = np.zeros((nofDataCols))
			means = np.zeros((nofDataCols))

			dataColi = -1
			for coli in dataCols:
				dataColi += 1
				print "Column", coli,
				means[dataColi] = np.mean(cols[coli])
				variances[dataColi] = np.var(cols[coli])
				stds[dataColi] = np.std(cols[coli])
				print "mean stdev skew kurt ", means[dataColi], stds[dataColi], scipy.stats.skew(cols[coli]), scipy.stats.kurtosis(cols[coli])

			# Autocorrelation functions on the entire series
			maxAcor = 0
			for coli in dataCols:
				if args.acf:
					# Autocorrelation time estimate (2009 Grossfield)
					# Compute correlation function for M values or until it reaches 0
					n = cols[coli].size # Sample size
					#M = args.Ms[coli] # Calculate up to point M
					M = 0 # Calculate up to point M
					if M == 0: M = n
		
					# Autocorrelation functions calculated with different methods
					corrFull, cut = CestGrossfield(M, cols[coli]) # Beyond the cut point it becomes noisy
					
					# Integrated autocorrelation time
					IAcFull = np.sum(corrFull)
					print "Integrated autocorrelation time full", IAcFull
					
					if maxAcor < IAcFull: maxAcor = int(IAcFull)
			print "Max full autocorr time", maxAcor

			# Get moving averages
			maCols = np.zeros((nofDataCols, (cols[0]).size ))
			colsMeanLines = np.ones((nofDataCols, (cols[0]).size ))
			dataColi = -1
			for coli in dataCols:
				dataColi += 1
				colsMeanLines[dataColi] = means[dataColi]
				if is_odd(maxAcor):
					maCols[dataColi] = moving_average(cols[coli], maxAcor)
				else:
					maCols[dataColi] = moving_average(cols[coli], maxAcor + 1)

			# Get moving average vs mean intersection
			firstIntersect = 0
			diff = prevDiff = 0.0
			dataColi = -1
			for coli in dataCols:
				dataColi += 1
				for i in range((cols[0]).size):
					if not np.isnan(maCols[dataColi][i]):
						diff = maCols[dataColi][i] - colsMeanLines[dataColi][i]
						if (diff * prevDiff) < 0:
							firstIntersect = i
							break
						prevDiff = diff

			# Get equilibration point
			eqPoint = 2 * firstIntersect
			print "Equilibration point", eqPoint

			# Get production run (Trim the series)
			trimCols = np.zeros((ncols, (cols[0][eqPoint:]).size ))
			for coli in range(ncols):
				trimCols[coli] = cols[coli][eqPoint:]


			# Get autocorrelation funtions on production runs
			dataColi = -1
			for coli in dataCols:
				dataColi += 1
				if args.acf:
					# Autocorrelation time estimate (2009 Grossfield)
					# Compute correlation function for M values or until it reaches 0
					n = trimCols[coli].size # Sample size
					#M = args.Ms[coli] # Calculate up to point M
					M = 0
					if M == 0: M = n
		
					# Autocorrelation functions calculated with different methods
					corr = np.zeros((7, M)) # (# of autocor functions, X size)
					corr[0], cut = CestGrossfield(M, trimCols[coli]) # Beyond the cut point it becomes noisy
					
					# Integrated autocorrelation time
					Iac = np.sum(corr[0])
					ESS = n / Iac
					print "Integrated autocorrelation time", Iac 
					print "Grossfield independent samples", ESS
		
					# Additional methods
					add_funcs = [autocorr1, autocorr2, autocorr3, autocorr4, autocorr5]
					for i in range(2, args.nofAddMethods):
						corr[i] = add_funcs[i](trimCols[coli], lags)
		
					# Standard error based on autocorrelation time fitting
					SE = stds[dataColi] * np.sqrt(Iac / float(n))
					print "Standard error", SE
		
					# Plots
					if(('acf' in args.makeplots) or ('all' in args.makeplots)):
						fignum += 1
						fig = plt.figure(fignum)
						fig.suptitle("Autocorrelation time function") 
						ax = plt.subplot(nofDataCols, 1, dataColi + 1)
						line, = ax.plot(corr[1], label='Col ' +  str(coli) + ' ' + ' f1', color='black')

						line.set_dashes([2,2,10,2])
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
					fig.suptitle(FNlist[li])
					ax = plt.subplot(nofDataCols, 1, dataColi+1)

					# Initial data
					ax.plot(cols[coli], label='Col ' +  str(coli) + ' real', color='gray')
					# Moving average
					ax.plot(maCols[dataColi], label='Col ' +  str(coli) + ' real', color='green')
					ax.plot(colsMeanLines[dataColi], label='Col ' +  str(coli) + ' real', color='cyan')
					# Trimmed data
					ax.plot( np.concatenate((np.zeros((eqPoint)), trimCols[coli])), label='Col ' +  str(coli) + ' real', color=colors[li])
	
					ax.plot(moving_average((1/avgAcc) * means[dataColi] *  cols[-1], maxAcor), color='blue')

					ax.set_xlabel(r'$\mathrm{t}$', fontsize=8)
					ax.set_ylabel(r'$\mathrm{f(t)}$', fontsize=8)

					if args.savefig:
						  figFN = 'temp.data.pdf'
						  plt.savefig(figFN, dpi=600, format='pdf')
	
if args.makeplots:
	plt.show()
	
	
	
	
	
	
	
	
	
	
	
	
