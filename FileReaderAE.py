import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm
from SampleDimensions import l_Q, w_Q, t_Q

l = np.delete(l_Q, 3)
w = np.delete(w_Q, 3)
t = np.delete(t_Q, 3)
a = w * t

AE_Q_Directory = r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\Ae\\Quasi"
dfAEQ = []

for file in os.scandir(AE_Q_Directory):
    dfAEQ.append(pd.read_csv(file))

for n, df in enumerate(dfAEQ):
    print(df.head())
    load = np.array(df["Load [N]"])
    E = np.array(df["E [eu]"])
    plt.plot(load, E, label=f'{n+1}')
    plt.show()