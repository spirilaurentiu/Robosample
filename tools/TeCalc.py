from statistics import mode
import numpy as np
import scipy.special

nof_replicas = 14
baseTemperature = 300.0
baseTdiff = 30.0

# temperatures = np.zeros(nof_replicas, dtype=np.float64)
# Tratio = (baseTemperature + baseTdiff) / baseTemperature
# for replIx in range(nof_replicas):
#     temperatures[replIx] = baseTemperature * (Tratio**replIx)
# print(temperatures)

# Temperatures distribution of REX
def REX_Ts_unkn(nof_replicas, base_temp, base_tdiff, dratio=0.0):
    """
    Calculate the temperatures for a replica exchange.
    """
    # Create the sequence of ratios: [Tratio, Tratio + dratio, Tratio + 2*dratio...]
    initial_ratio = (base_temp + base_tdiff) / base_temp
    ratios = initial_ratio + np.arange(nof_replicas - 1) * dratio
    
    # Prepend 1.0 for the first replica (which stays at base_temp)
    full_ratios = np.concatenate(([1.0], ratios))
    
    return base_temp * np.cumprod(full_ratios)

# temperatures = REX_Ts_unkn(nof_replicas, baseTemperature, baseTdiff, dratio=0.08)
# for replIx in range(nof_replicas):
#     print(f"Replica {replIx}: {temperatures[replIx]:.2f} K")

# Quadratic distribution of temperatures for REX
def REX_Ts_test(nof_replicas, base_temp, base_tdiff):
    """Calculate Ts for replica exchange using a quadratic."""

    T = np.zeros(nof_replicas, dtype=np.float64)

    replIx = 0
    for dix in range(0, nof_replicas, 1):
        T[replIx] = (7*(dix**2.5)) + (0*dix) + base_temp
        replIx += 1

    return T
#

def printTabel(headers, data):
    """
    Docstring for printTabel
    
    :param headers: list of strings representing the column headers
    :param data: list of arrays representing columns
    """
    # Print headers
    print("\t".join(headers))
    
    min_rows = min(len(col) for col in data)
    max_rows = max(len(col) for col in data)
    
    # Print data rows
    for rowIx in range(max_rows):
        row = []
        for col in data:
            if rowIx < len(col):
                row.append(col[rowIx])
            else:
                row.append("None")  # Placeholder for missing data

        for value in row:
            print(value, end="\t")
        print()  # Newline after each row

def REX_Ts_exp(n_replicas, T_min, T_max):
    """
    Calculates temperatures using a geometric progression.
    The ratio between adjacent temperatures is constant.
    """
    # np.geomspace creates a geometric progression from T_min to T_max
    Ts = np.geomspace(T_min, T_max, num=n_replicas)
    
    # To use your printTabel, we need to calculate the "step" progress
    steps = np.arange(n_replicas)
    
    headers = ["Step", "Ts"]
    data = [steps, Ts]
    #printTabel(headers, data)
    
    return Ts

def REX_Ts_logit(n_replicas, T_min, T_max):

    x = np.linspace(0.05, 0.95, n_replicas)
    # x[0] = 0.000001  # Avoid logit(0)
    # x[-1] = 0.999999  # Avoid logit(1)
    logitF = scipy.special.logit(x)

    # Normalize it so it starts at 0 and ends at 1
    logitnorm = logitF - logitF[0]  # Shift to start at 0
    logitnorm = logitnorm / logitnorm[-1]  # Scale to end at 1

    Ts = T_min + (T_max - T_min) * logitnorm
    headers = ["x", "logit", "logit_norm", "Ts"]
    data = [x, logitF, logitnorm, Ts]
    #printTabel(headers, data)

    return Ts
#

def REX_Ts_probit(n_replicas, T_min, T_max):
    # Probit is the inverse of the Normal CDF (ndtri)
    # We use 0.05 to 0.95 to avoid -inf and +inf
    x = np.linspace(0.05, 0.95, n_replicas)
    probitF = scipy.special.ndtri(x)

    # Normalize 0 to 1
    norm = (probitF - probitF[0]) / (probitF[-1] - probitF[0])
    Ts = T_min + (T_max - T_min) * norm
    
    printTabel(["x", "probit", "norm", "Ts"], [x, probitF, norm, Ts])
    return Ts

def REX_Ts_arcsin(n_replicas, T_min, T_max):
    # Arcsin maps -1 to 1, we map our x to that range
    x = np.linspace(-1, 1, n_replicas)
    arcsinF = np.arcsin(x)

    # Normalize 0 to 1
    norm = (arcsinF - arcsinF[0]) / (arcsinF[-1] - arcsinF[0])
    Ts = T_min + (T_max - T_min) * norm
    
    printTabel(["x", "arcsin", "norm", "Ts"], [x, arcsinF, norm, Ts])
    return Ts

temperatures1 = REX_Ts_arcsin(7, 300, 500)
temperatures2 = REX_Ts_exp(8, 500, 1000)

# temperatures = np.concatenate((temperatures1, temperatures2[1:]))  # Combine, skipping the first of the second set to avoid duplication
# print("\nFinal combined temperatures:")
# for replIx in range(len(temperatures)):
#     print(f"Replica {replIx}: {temperatures[replIx]:.2f} K")

def REX_Ts(mode, n_replicas, T_min, T_max, **kwargs):
    """
    Master function to generate temperature distributions.
    Modes: 'logit', 'probit', 'arcsin', 'exp', 'power'
    """
    steps = np.arange(n_replicas)

    if mode == 'gompertz':
        # b controls displacement, c controls growth rate
        b, c = kwargs.get('b', 5), kwargs.get('c', 10)
        raw_vals = np.exp(-b * np.exp(-c * steps))
        
    elif mode == 'sine':
        # Using the first quarter of a sine wave (0 to 90 degrees)
        raw_vals = 1 - np.cos(steps * (np.pi / 2))
        
    elif mode == 'asym_sig':
        # k controls the "steepness"
        # shift controls where the midpoint is (0.5 is symmetric)
        k = kwargs.get('k', 10)
        shift = kwargs.get('shift', 0.7) # 0.7 makes it stay slow longer
        raw_vals = 1 / (1 + np.exp(-k * (steps - shift)))
            
    elif mode == 'exp':
        baseDiff = kwargs.get('baseDiff', 10.0)
        ratio = (baseTemperature + baseDiff) / baseTemperature
        raw_vals = np.zeros(n_replicas, dtype=np.float64)
        for replIx in range(n_replicas):
            raw_vals[replIx] = baseTemperature * (ratio ** replIx)
        return raw_vals
        
    elif mode in ['logit', 'probit', 'arcsin']:
        # Define the internal mapping range to avoid +/- infinity
        # Logit/Probit need (0, 1), Arcsin needs (-1, 1)
        if mode == 'arcsin':
            x = np.linspace(-0.95, 0.95, n_replicas)
            raw_vals = np.arcsin(x)
        elif mode == 'probit':
            x = np.linspace(0.05, 0.95, n_replicas)
            raw_vals = scipy.special.ndtri(x)
        else: # logit
            x = np.linspace(0.05, 0.95, n_replicas)
            raw_vals = scipy.special.logit(x)

    elif mode == 'power':
        p = kwargs.get('power', 2.0)
        x = np.linspace(0, 1, n_replicas)
        raw_vals = x ** p
        Ts = T_min + (T_max - T_min) * raw_vals

    else:
        raise ValueError(f"Unknown mode: {mode}")
                
    # Normalize and Scale
    norm = (raw_vals - raw_vals[0]) / (raw_vals[-1] - raw_vals[0])
    Ts = T_min + (T_max - T_min) * norm

    # Use your custom print function
    headers = ["Index", "Raw_Val", "Temp"]
    data = [steps, raw_vals, Ts]
    printTabel(headers, data)
    
    return Ts
#

temperatures_exp = REX_Ts('exp', 14, 300, 1000, baseDiff=80)
teperatures_pow = REX_Ts('power', 14, 300, 1000, power=2.5)
temperatures_logit = REX_Ts('logit', 14, 300, 1000)
temperatures_probit = REX_Ts('probit', 14, 300, 1000)
temperatures_arcsin = REX_Ts('arcsin', 14, 300, 1000)
temperatures_gompertz = REX_Ts('gompertz', 14, 300, 1000)
temperatures_sin = REX_Ts('sine', 14, 300, 1000)
temperatures_asym_sig = REX_Ts('asym_sig', 14, 300, 1000, k=10, shift=0.7)
printTabel(["Exp", "Power", "Logit", "Probit", "Arcsin", "Gompertz", "Sine", "Asym_Sig"],
           [temperatures_exp, teperatures_pow, temperatures_logit, temperatures_probit, temperatures_arcsin, temperatures_gompertz, temperatures_sin, temperatures_asym_sig])
