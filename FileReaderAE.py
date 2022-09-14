import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_Q, w_Q, t_Q

l = np.delete(l_Q, 3)
w = np.delete(w_Q, 3)
t = np.delete(t_Q, 3)
a = w * t

AE_Q_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\Ae\\Quasi"
dfAEQ = []
peaks = [76, 77, 98, 98, 82, 66, 68, 84, 74]
initiationloads = []
failureloads = []


for file in os.scandir(AE_Q_Directory):
    dfAEQ.append(pd.read_csv(file))

def closest(lst, K):
    lst = np.asarray(lst)
    idx = (np.abs(lst - K)).argmin()
    return idx

for n, df in enumerate(dfAEQ):

    t = list(df["Time [s]"])
    load = list(df["Load [N]"])
    E = np.array(df["E [eu]"])
    cumeg = [E[0]]
    for j, E_point in enumerate(E[1:]):
        new_cumeg = cumeg[j] + E_point
        cumeg.append(new_cumeg)
    print(len(cumeg), len(E))

    t_min = closest(t, 60)
    t_index = closest(t, peaks[n])
    f_max = max(load)/w[n]
    forceatpeak = load[t_index]/w[n]

    initiationloads.append(forceatpeak)
    failureloads.append(f_max)

    plt.plot(t[t_min:], E[t_min:], 'g', label=f'Sample {n+1} - # of hits')
    #plt.plot(t[t_min:], cumeg[t_min:], 'y', label=f'Sample {n + 1} - c. energy')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [eu]')
    plt.axvline(t[t_index])
    plt.grid()
    plt.legend()
    plt.show()


print(initiationloads)
print(failureloads)
mu, std = norm.fit(initiationloads)
mu2, std2 = norm.fit(failureloads)
print(mu, mu2)

