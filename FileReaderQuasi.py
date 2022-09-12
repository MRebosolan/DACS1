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


MTS_UD_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\Quasi"
dfMTSQ = []
UTSlist = []


for file in os.scandir(MTS_UD_Directory):
    dfMTSQ.append(pd.read_csv(file))

for n, df in enumerate(dfMTSQ):
    sigmaxx = np.array(df["Load [N]"]/a[n])
    exx = np.array(df["Displacement [mm]"]/l[n])
    UTS = max(sigmaxx)
    UTSlist.append(UTS)

print(UTSlist*t)



mu, std = norm.fit(UTSlist)
print(mu)
