import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_90, w_90, t_90


DIC_90_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DIC\\90"
dfDIC90 = []
E2values = []
Ytvalues = []
v21values = []
G12values = []
l = l_90[2:]
w = w_90[2:]
t = t_90[2:]
a = w * t

for file in os.scandir(DIC_90_Directory):
    dfDIC90.append(pd.read_csv(file, sep=";", skiprows=1))


for n, df in enumerate(dfDIC90):
    df = df.dropna()
    exx = np.array(df["exx [1] - engr."])
    eyy = np.array(df["eyy [1] - engr."])
    sigmaxx = np.array(df["Force"]/a[n])
    G12 = np.array(sigmaxx[30:50] / (2 * exx[30:50] * eyy[30:50]))
    v21 = - eyy[20:50]/exx[20:50]
    for value in v21:
        v21values.append(abs(value))
    for value2 in G12:
        G12values.append(value2)

mu, std = norm.fit(v21values)
mu2, std2 = norm.fit(G12values)
print(mu2)


