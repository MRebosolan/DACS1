import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_UD, w_UD, t_UD
from reliability.Fitters import Fit_Weibull_2P, Fit_Weibull_3P
from reliability.Probability_plotting import plot_points


l = np.delete(l_UD, [6, 9])
w = np.delete(w_UD, [6, 9])
t = np.delete(t_UD, [6, 9])
E1values = []
v12values = []
G12values = []
a = w * t

print(a)




DIC_UD_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DIC\\UD"
dfDICUD = []

for file in os.scandir(DIC_UD_Directory):
    dfDICUD.append(pd.read_csv(file, sep=";", skiprows=1))

for n, df in enumerate(dfDICUD):
    if n != 8:
        df = df.dropna()
        exx = np.array(df["exx [1] - engr."])
        eyy = np.array(df["eyy [1] - engr."])

        sigmaxx = np.array(df["Force"]/(a[n]))
        G12 = np.array(sigmaxx[30:50]/(2*(exx[30:50] - eyy[30:50])))

        E1 = sigmaxx[20:50]/exx[20:50]
        v12 = -eyy[20:50]/exx[20:50]
        #plt.plot(exx, sigmaxx,  label=f'{n+1}')
        for value in v12:
            v12values.append(value)
        for value2 in E1:
            E1values.append(value2)
        for value3 in G12:
            G12values.append(G12)


mu, std = norm.fit(v12values)
mu2, std2 = norm.fit(E1values)
mu3, std3 = norm.fit(G12values)

data = E1values
weibull_fit = Fit_Weibull_3P(failures=data,show_probability_plot=False,print_results=False)
weibull_fit.gamma = 120
weibull_fit.distribution.PDF(label='Fitted Distribution',color='steelblue')
#plot_points(failures=data,func='PDF',label='failure data',color='red',alpha=0.7)

plt.legend()
plt.show()




