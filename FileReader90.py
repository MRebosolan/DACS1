import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_90, w_90, t_90
from reliability.Fitters import Fit_Weibull_2P


MTS_90_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\90"
dfMTS90 = []
E2values = []
Ytvalues = []
l = l_90[2:]
w = w_90[2:]
t = t_90[2:]
a = w * t



for file in os.scandir(MTS_90_Directory):
    dfMTS90.append(pd.read_csv(file))


for n, df in enumerate(dfMTS90):

    df["eyy"] = df["Displacement [mm]"]/l[n]
    df["sigmayy"] = df["Load [N]"]/a[n]
    eyy = np.array(df["eyy"])
    sigmayy = np.array(df["sigmayy"])
    Ytvalues.append(max(sigmayy))
    E2 = sigmayy[1000:5000]/eyy[1000:5000]
    for value in E2:
        E2values.append(value)

mu, std = norm.fit(E2values)
mu2, std2 = norm.fit(Ytvalues)
print(mu2, std2)

data = Ytvalues
weibull_fit = Fit_Weibull_2P(failures=data,show_probability_plot=False,print_results=False)
weibull_fit.distribution.PDF(label='Fitted Distribution - Yt',color='steelblue')

plt.xlabel("Yt (MPa)")
plt.xlim((60, 140))
plt.grid()
plt.legend()
plt.show()


