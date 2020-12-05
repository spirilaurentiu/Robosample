import numpy as np

def func0(series):
	result = np.zeros((series.size))
	for i in range(series.size):
		result[i] = (np.mean(series[:i]) / np.mean(series[i:]))
	return result
#

def func1(series):
	result = np.zeros((series.size))
	for i in range(series.size):
		result[i] = (np.mean(series[:i]) / np.mean(series))
	return result
#

def func2(series):
	result = np.zeros((series.size))
	for i in range(series.size):
		result[i] = np.log(np.sum(np.exp(series[:i])))
	return result
#

