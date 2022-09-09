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
l = l_90[2:]
w = w_90[2:]
t = t_90[2:]
a = w * t

for file in os.scandir(DIC_90_Directory):
    dfDIC90.append(pd.read_csv(file, sep=";", skiprows=1))


for df in dfDIC90:
    df = df.dropna()
    exx = np.array(df["exx [1] - engr."])
    eyy = np.array(df["eyy [1] - engr."])
    v21 = - eyy[20:50]/exx[20:50]
    for value in v21:
        v21values.append(abs(value))

mu, std = norm.fit(v21values)
print(v21values)
print(mu)


