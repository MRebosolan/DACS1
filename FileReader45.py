import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_45, w_45, t_45

MTS_45_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\45"
dfMTS45 = []
S12values = []

l = np.delete(l_45, [0, 3, 4, 8])
w = np.delete(w_45, [0, 3, 4, 8])
t = np.delete(t_45, [0, 3, 4, 8])
a = w * t


for file in os.scandir(MTS_45_Directory):
    dfMTS45.append(pd.read_csv(file))

for n, df in enumerate(dfMTS45):
    load = np.array(df["Load [N]"])
    exx = np.array(df["Displacement [mm]"])/l[n]
    sigmaxx = np.array(df["Load [N]"])/a[n]
    S12 = max(sigmaxx)
    S12values.append(S12)

plt.show()

print(S12values)