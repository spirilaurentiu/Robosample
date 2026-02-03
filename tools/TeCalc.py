import numpy as np

nof_replicas = 21
baseTemperature = 300.0
baseTdiff = 20.0

temperatures = np.zeros(nof_replicas, dtype=np.float64)
Tratio = (baseTemperature + baseTdiff) / baseTemperature

for replIx in range(nof_replicas):
    temperatures[replIx] = baseTemperature * (Tratio**replIx)
print(temperatures)
