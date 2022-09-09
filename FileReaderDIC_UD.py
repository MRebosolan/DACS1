import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_UD, w_UD, t_UD




l = np.delete(l_UD, [6, 10])
w = np.delete(w_UD, [6, 10])
t = np.delete(t_UD, [6, 10])
E1values = []
v12values = []
G12values = []
a = w * t




DIC_UD_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DIC\\UD"
dfDICUD = []

for file in os.scandir(DIC_UD_Directory):
    dfDICUD.append(pd.read_csv(file, sep=";", skiprows=1))

for n, df in enumerate(dfDICUD):

    df = df.dropna()
    exx = np.array(df["exx [1] - engr."])
    eyy = np.array(df["eyy [1] - engr."])

    sigmaxx = np.array(df["Force"]/(a[n]))
    G12 = np.array(sigmaxx[30:50]/(2*(exx[30:50] - eyy[30:50])))
    E1 = sigmaxx[30:50]/exx[30:50]
    v12 = -eyy[30:50]/exx[30:50]
    for value in v12:
        v12values.append(value)
    for value2 in E1:
        E1values.append(value2)
    for value3 in G12:
        G12values.append(G12)


mu, std = norm.fit(v12values)
mu2, std2 = norm.fit(E1values)
mu3, std3 = norm.fit(G12values)
print(mu3)




