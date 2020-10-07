import sys, os, glob
import numpy as np
import scipy
import scipy.stats

def moving_average(series, window_size):
        if np.size(series.shape) > 1:
                print("Moving average function error: Implemented only for 1D data.")
                return 1
        if (window_size % 2 != 1):
		window_size = window_size + 1
                print("Moving average function warning: Window size increse by 1.")

        shift = (window_size - 1) / 2
        #retSeries = np.ones((series.shape[0] - (2*shift)))
        retSeries = np.zeros((series.shape[0])) * np.nan
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

def moving_difference(series, window_size):
        if np.size(series.shape) > 1:
                print("Moving difference function error: Implemented only for 1D data.")
                return 1
        if (window_size % 2 != 1):
                print("Moving difference function error: Implemented only for odd window sizes.")
                return 1
        shift = (window_size - 1) / 2
        retSeries = np.zeros((series.shape[0]))
        for i in range(shift, (retSeries.shape[0] - shift)):
                retSeries[i] = (series[i + shift] - series[i - shift]) / window_size
        return retSeries
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

        return resultF, cut
#

# Functions from https://stackoverflow.com/users/2005415/jason
def autocorr1(argM, x):
    '''numpy.corrcoef, partial'''

    lags = np.arange(0, argM)
    corr=[1. if l==0 else numpy.corrcoef(x[l:],x[:-l])[0][1] for l in lags]
    return numpy.array(corr)

def autocorr2(argM, x):
    '''manualy compute, non partial'''

    lags = np.arange(0, argM)
    mean=numpy.mean(x)
    var=numpy.var(x)
    xp=x-mean
    corr=[1. if l==0 else numpy.sum(xp[l:]*xp[:-l])/len(x)/var for l in lags]

    return numpy.array(corr)

def autocorr3(argM, x):
    '''fft, pad 0s, non partial'''

    lags = np.arange(0, argM)
    n=len(x)
    # pad 0s to 2n-1
    ext_size=2*n-1
    # nearest power of 2
    fsize=2**numpy.ceil(numpy.log2(ext_size)).astype('int')

    xp=x-numpy.mean(x)
    var=numpy.var(x)

    # do fft and ifft
    cf=numpy.fft.fft(xp,fsize)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real
    corr=corr/var/n

    return corr[:len(lags)]

def autocorr4(argM, x):
    '''fft, don't pad 0s, non partial'''
    lags = np.arange(0, argM)
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean

    cf=numpy.fft.fft(xp)
    sf=cf.conjugate()*cf
    corr=numpy.fft.ifft(sf).real/var/len(x)

    return corr[:len(lags)]

def autocorr5(argM, x):
    '''numpy.correlate, non partial'''
    lags = np.arange(0, argM)
    mean=x.mean()
    var=numpy.var(x)
    xp=x-mean
    corr=numpy.correlate(xp,xp,'full')[len(x)-1:]/var/len(x)

    return corr[:len(lags)]
#
