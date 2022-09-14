import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_45, w_45, t_45
from reliability.Fitters import Fit_Weibull_2P


DIC_45_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DIC\\45"
dfDIC45 = []


l = np.delete(l_45, [0, 3, 4, 8])
w = np.delete(w_45, [0, 3, 4, 8])
t = np.delete(t_45, [0, 3, 4, 8])
a = w * t

G12values = []

for file in os.scandir(DIC_45_Directory):
    dfDIC45.append(pd.read_csv(file, sep=";", skiprows=1))

for n, df in enumerate(dfDIC45):

    exx = np.array(df["exx [1] - engr."])
    eyy = np.array(df["eyy [1] - engr."])
    sigmaxx = np.array(df["Force"]) / a[n]
    G12 = sigmaxx[30:50]/(2*(exx[30:50]-eyy[30:50]))
    for value in G12:
        G12values.append(value)

mu, std = norm.fit(G12values)
print(mu, std)


data = G12values
weibull_fit = Fit_Weibull_2P(failures=data, show_probability_plot=False, print_results=False)
weibull_fit.distribution.PDF(label='Fitted Distribution - G12', color='steelblue')

plt.xlabel("G12 (GPa)")
plt.xlim((1, 6))
plt.grid()
plt.legend()
plt.show()
