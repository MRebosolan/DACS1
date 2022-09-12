import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_Q, w_Q, t_Q

MTS_Q_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\45"
dfMTSQ = []
failureloads = []

l = np.delete(l_Q, 3)
w = np.delete(w_Q, 3)
t = np.delete(t_Q, 3)
a = w * t

for file in os.scandir(MTS_Q_Directory):
    dfMTSQ.append(pd.read_csv(file))

for df in dfMTSQ:
    print(df.head())