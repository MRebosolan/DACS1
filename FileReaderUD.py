import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm, exponweib
from SampleDimensions import l_UD, w_UD, t_UD
from reliability.Fitters import Fit_Weibull_2P, Fit_Weibull_3P



#APPROACH TO OBTAIN E1
#IT CAN BE OBAINED FROM DIC_UD AND MTS_UD
#AVERAGE t,w,l FOR EACH SAMPLE
#EXTRACT STRAINS AND STRESSES FROM LOADS AND DISPLACEMENTS (FOR MTS)





# l = np.array([250.17,250.16,250.15,250.14,250.15,250.16,250.15,250.13,250.11,250.16,250.17])
# a = np.array([2,2,2,2,2,2,2,2,2,2,2]) * 0.00001

l = np.delete(l_UD, [6, 9])
w = np.delete(w_UD, [6, 9])
t = np.delete(t_UD, [6, 9])
a = w * t




MTS_UD_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\UD"
dfMTSUD = []
E1values = []
Xtvalues = []


for file in os.scandir(MTS_UD_Directory):
    dfMTSUD.append(pd.read_csv(file))

for i, df in enumerate(dfMTSUD):
    df["exx"] = df["Displacement [mm]"]/l[i]
    df["sigmaxx"] = df["Load [N]"]/(a[i])
    exx = np.array(df["exx"])
    sigmaxx = np.array(df["sigmaxx"])

    #plt.plot(exx, sigmaxx, label=f"{i+1}")

    Xt = max(sigmaxx)
    Xtvalues.append(Xt)
    E1 = sigmaxx[1000:5000]/exx[1000:5000]
    for value in E1:
        E1values.append(value)



# mu, std = norm.fit(E1values)
# mu2, std2 = norm.fit(Xtvalues)
x = np.linspace(min(E1values), max(E1values), 100000)
# x2 = np.linspace(min(Xtvalues), max(Xtvalues), 100000)
# p = norm.pdf(x, mu, std)
# p2 = norm.pdf(x2, mu2, std2)
# plt.legend()
# #plt.plot(x, p)
# #plt.plot(x2, p2)
# plt.show()
# print(mu2, std2)

# exp1, k1, loc1, lam1 = exponweib.fit(E1values, f0=1)
# weibullpdf = exponweib.pdf(x, exp1, k1, loc1, lam1)
# print(exp1, k1, loc1, lam1)


data = Xtvalues
weibull_fit = Fit_Weibull_2P(failures=data,show_probability_plot=False,print_results=False)
weibull_fit.distribution.PDF(label='Fitted Distribution - Xt',color='steelblue')

plt.xlabel("Xt (MPa)")
plt.xlim((1250, 2250))
plt.grid()
plt.legend()
plt.show()





