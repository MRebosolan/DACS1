import numpy as np
from math import *

E1 = 140.03E3
E2 = 7.72E3
v12 = 0.35
v21 =v12 * (E2/E1)
G12 = 4.685E3
S12 = 79
S23 = 20
Xt = 1950
Xc = 1480
Yt = 107
Yc = 220
t = 0.125

beta = 1.04E-8
G1c = 0.258
G2c = 1.080


phi = (48*G2c)/(pi*t)
phi2 = (24*G2c)/(pi*t)

S12insitu = sqrt((sqrt(1 + beta*phi*G12**2)-1)/(3*beta*G12))
S12insitu2 = sqrt((sqrt(1 + beta*phi2*G12**2)-1)/(3*beta*G12))
delta = 2*((1/E2) - (v21**2 / E1))
Ytinsitu = sqrt((8*G1c)/(pi*t*delta))
Ytinsitu2 = 1.79*sqrt(G1c/(pi*t*delta))

print(S12insitu, S12insitu2)