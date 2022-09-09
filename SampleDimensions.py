import pandas as pd
import numpy as np
from pathlib import Path


fileUD = Path(r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DimensionsUD.csv")
file90 = Path(r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\Dimensions90.csv")
file45 = Path(r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\Dimensions45.csv")
fileQ = Path(r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DimensionsQuasi.csv")



dfUD = pd.read_csv(fileUD)
df90 = pd.read_csv(file90)
df45 = pd.read_csv(file45)
dfQ = pd.read_csv(fileQ)


l_UD = np.array(dfUD["L Left"] + dfUD["L Right"])/2
w_UD = np.array(dfUD["W Top"] + dfUD["W Mid"] + dfUD["W Bottom"])/3
t_UD = np.array(dfUD["T Top-left"] + dfUD["T Top-right"] + dfUD["T Mid-left"] + dfUD["T Mid-right"] + dfUD["T Bottom-left"] + dfUD["T Bottom-right"])/6

l_90 = np.array(df90["L Left"] + df90["L Right"])/2
w_90 = np.array(df90["W Top"] + df90["W Mid"] + df90["W Bottom"])/3
t_90 = np.array(df90["T Top-left"] + df90["T Top-right"] + df90["T Mid-left"] + df90["T Mid-right"] + df90["T Bottom-left"] + df90["T Bottom-right"])/6

l_45 = np.array(df45["L Left"] + df45["L Right"])/2
w_45 = np.array(df45["W Top"] + df45["W Mid"] + df45["W Bottom"])/3
t_45 = np.array(df45["T Top-left"] + df45["T Top-right"] + df45["T Mid-left"] + df45["T Mid-right"] + df45["T Bottom-left"] + df45["T Bottom-right"])/6

l_Q = np.array(dfQ["L Left"] + dfQ["L Right"])/2
w_Q = np.array(dfQ["W Top"] + dfQ["W Mid"] + dfQ["W Bottom"])/3
t_Q = np.array(dfQ["T Top-left"] + dfQ["T Top-right"] + dfQ["T Mid-left"] + dfQ["T Mid-right"] + dfQ["T Bottom-left"] + dfQ["T Bottom-right"])/6



