import pandas as pd
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import os


#APPROACH TO OBTAIN E1
#IT CAN BE OBAINED FROM DIC_UD AND MTS_UD
#AVERAGE t,w,l FOR EACH SAMPLE
#EXTRACT STRAINS AND STRESSES FROM LOADS AND DISPLACEMENTS (FOR MTS)





l = np.array([250.17,250.16,250.15,250.14,250.15,250.16,250.15,250.13,250.11,250.16,250.17])
a = np.array([2,2,2,2,2,2,2,2,2,2,2]) * 0.00001



MTS_UD_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\MTS\\UD"
dfMTSUD = []
E1values = []

for file in os.scandir(MTS_UD_Directory):
    dfMTSUD.append(pd.read_csv(file))

for i, df in enumerate(dfMTSUD):
    df["exx"] = df["Displacement [mm]"]/l[i]
    df["sigmaxx"] = df["Load [N]"]/(a[i])
    exx = np.array(df["exx"])
    sigmaxx = np.array(df["sigmaxx"])
    plt.xlim((0, 0.014))
    plt.ylim(-100, 1.5E9)
    plt.plot(exx, sigmaxx)
    E1 = sigmaxx[10000]/exx[10000]
    print(E1)
    plt.show()

print(E1values)


