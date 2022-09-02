import pandas as pd
import numpy as np
from pathlib import Path


file = Path(r"C:\\Users\\nxf92804\\Desktop\\AE4ASM109 Data\\AssignmentCode\\ASM109_2021_data\\DimensionsUD.csv")

df = pd.read_csv(file)
print(df)