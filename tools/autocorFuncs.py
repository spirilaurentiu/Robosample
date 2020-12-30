import sys, os, glob
import numpy as np
import scipy
import scipy.stats

def is_odd(num):
	return num % 2 != 0
#

def printOneLiner1D(series):
	''' Print Numpy array on one line'''
	for i in range(series.shape[0]):
		print (series[i]),
	print()
#

def print2D(M):
	''' Print Numpy 2D array on row per line'''
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			print (M[i, j]),
		print()
#

def print2DAndExit(M):
	for i in range(M.shape[0]):
		for j in range(M.shape[1]):
			print (M[i, j]),
		print()
	exit(0)
#

def moving_average(series, window_size):
	if np.size(series.shape) > 1:
		print("Moving average function error: Implemented only for 1D data.")
		return 1
	if (window_size % 2 != 1):
		window_size = window_size + 1
		print("Moving average function warning: Window size increse by 1.")

	shift = int((window_size - 1) / 2)
	#retSeries = np.ones((series.shape[0] - (2*shift)))
	retSeries = np.zeros((series.shape[0])) 
	for i in range(shift, (retSeries.shape[0] - shift)):
		retSeries[i] = np.mean(series[i - shift: i + shift + 1])
	return retSeries
#

def cumulative_average(series):
	if np.size(series.shape) > 1:
		print("Cumulative average function error: Implemented only for 1D data.")
		return 1
	retSeries = np.zeros((series.shape[0]))
	retSeries[0] = series[0]
	for i in range(1, retSeries.shape[0]):
		retSeries[i] = np.mean(series[0:i])
	return retSeries
#

def difference_quotient(series, window_size):
	if np.size(series.shape) > 1:
		print("Difference quotient function error: Implemented only for 1D data.")
		return 1
	if (window_size % 2 != 1):
		window_size = window_size + 1
		print("Difference quotient function warning: Window size increse by 1.")
	shift = int((window_size - 1) / 2)
	retSeries = np.zeros((series.shape[0]))
	for i in range(shift, (retSeries.shape[0] - shift)):
		retSeries[i] = (series[i + shift] - series[i - shift]) / window_size
	return retSeries
#

def moving_difference(series, window_size):
	if np.size(series.shape) > 1:
		print("Moving difference function error: Implemented only for 1D data.")
		return 1
	if (window_size % 2 != 1):
		window_size = window_size + 1
		print("Moving difference function warning: Window size increse by 1.")
	shift = int((window_size - 1) / 2)
	retSeries = np.zeros((series.shape[0]))
	for i in range(shift, (retSeries.shape[0] - shift)):
		retSeries[i] = (series[i + shift] - series[i - shift]) 
	return retSeries
#

def intersections(series1, series2):
	'''Intersection of two series'''
	minSize = np.min([series1.size, series2.size])
	intersections = []
	diff = prevDiff = 0.0
	for i in range(minSize):
		if not np.isnan(series1[i]):
			diff = series1[i] - series2[i]
			if (diff * prevDiff) < 0:
				intersections.append(i)
				#break
			prevDiff = diff
	return np.array(intersections)
#


# Take the product (f(t))(f(t+lag)) and average over (n - lag) points
# Normalize. The variance is also the correlation at time 0 and
# should be highesta
def CestGrossfield(argM, series):
	# Define a zero
	tiny = 0.0000001

	# Moments
	n = series.size
	mean = np.mean(series)
	variance = np.var(series)

	# Vars
	lags = np.arange(0, argM)
	cut = 0
	resultF = np.zeros((argM))

	for lag in lags:

		resultF[lag] = np.mean( \
			(series[0 : (n - lag)] - mean) * \
			(series[lag : n]       - mean) ) \
			/ variance

		# If the correlation function reaches 0 the rest is considered noise
		if resultF[lag] < tiny:
			cut = lag
			break

	# Integrated autocorrelation time and effective sample size
	integratedAC = np.sum(resultF)
	effSampleSize = n / integratedAC	

	return resultF, cut, integratedAC, effSampleSize
#

autocorrLabels = ['FFT Partial', 'FFT Nonpartial', 'NumPy Correlate Nonpartial', 'NumPy Correlate Partial', 'Manual']

# Functions from https://stackoverflow.com/users/2005415/jason
def autocorrDetail(argM, series):
	'''fft, don't pad 0s, non partial'''

	# Define a zero
	tiny = 0.0000001

	n = series.size

	if argM == 0:
		argM = n
	lags = np.arange(0, argM)
	mean = series.mean()
	var = np.var(series)
	seriesp = series - mean

	cf = np.fft.fft(seriesp)
	freq = np.fft.fftfreq(n)
	#return (freq, cf)
	sf = cf.conjugate()*cf
	#return (freq, sf)
	invSf = np.fft.ifft(sf)

	resultF = invSf.real/var/len(series)

	return (freq, cf, sf, invSf, resultF)

	#cut = 0
	#for i in range(resultF.shape[0]):
	#	if resultF[i] < tiny:
	#		cut = i
	#		break
	#resultF[cut:] = 0


	#print(resultF)
	#return resultF[:len(lags)]

	# Integrated autocorrelation time and effective sample size
	#integratedAC = np.sum(resultF)
	#effSampleSize = n / integratedAC	

	#return resultF, cut, integratedAC, effSampleSize
#

# Functions from https://stackoverflow.com/users/2005415/jason
def autocorr1(argM, series):
	'''fft, don't pad 0s, non partial'''

	# Define a zero
	tiny = 0.0000001

	n = series.size

	lags = np.arange(0, argM)
	mean = series.mean()
	var = np.var(series)
	seriesp = series - mean

	cf = np.fft.fft(seriesp)
	sf = cf.conjugate()*cf
	resultF = np.fft.ifft(sf).real/var/len(series)

	cut = 0
	for i in range(resultF.shape[0]):
		if resultF[i] < tiny:
			cut = i
			break
	resultF[cut:] = 0


	#print(resultF)
	#return resultF[:len(lags)]

	# Integrated autocorrelation time and effective sample size
	integratedAC = np.sum(resultF)
	effSampleSize = n / integratedAC	

	return resultF, cut, integratedAC, effSampleSize
#

def autocorr2(argM, x):
	'''fft, pad 0s, non partial'''

	lags = np.arange(0, argM)
	n=len(x)
	# pad 0s to 2n-1
	ext_size=2*n-1
	# nearest power of 2
	fsize=2**np.ceil(np.log2(ext_size)).astype('int')

	xp=x-np.mean(x)
	var=np.var(x)

	# do fft and ifft
	cf=np.fft.fft(xp,fsize)
	sf=cf.conjugate()*cf
	corr=np.fft.ifft(sf).real
	corr=corr/var/n

	return corr[:len(lags)]

def autocorr3(argM, x):
	'''numpy.corrcoef, partial'''

	lags = np.arange(0, argM)
	corr=[1. if l==0 else np.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
	return np.array(corr)

def autocorr4(argM, x):
	'''numpy.correlate, non partial'''
	lags = np.arange(0, argM)
	mean=x.mean()
	var=np.var(x)
	xp=x-mean
	corr=np.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

	return corr[:len(lags)]
#

def autocorr5(argM, x):
	'''manualy compute, non partial'''

	lags = np.arange(0, argM)
	mean=np.mean(x)
	var=np.var(x)
	xp=x-mean
	corr=[1. if l==0 else np.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

	return np.array(corr)

