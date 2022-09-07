import numpy as np
from math import *
import matplotlib.pyplot as plt
from Functions import *

angles = []

n = 2

for i in range(n):
    angles.append(0)
    angles.append(90)
    angles.append(45)
    angles.append(-45)
    angles.append(-45)
    angles.append(45)
    angles.append(90)
    angles.append(0)

E1 = 140.03E3
E2 = 7.72E3
v12 = 0.28
G12 = 4.685E3
S12 = 79
Xt = 1950
Xc = 1200
Yt = 48
Yc = 35

t = 0.125

zlocations = np.arange(-t*len(angles)/2, t*(len(angles)+1)/2, t)

laminaarray = []

for n in range(len(angles)):
    ply = lamina(E1, E2, v12, G12, radians(angles[n]), zlocations[n], zlocations[n+1], Xt, Xc, Yt, Yc, S12)
    laminaarray.append(ply)

laminate1 = laminate(laminaarray)


Nxrange = np.arange(-1000, 1000, 25)
Nyrange = np.arange(0, 1000, 25)

envelopepoints = []
loadcombinations = {}  #KEYS = ENVELOPE COORDINATES, VALUES = FAILED TRUE/FALSE

#PLYSTRAINS --> CALCULATESTRESSES --> PRINCIPLASTRESSES2

for Nx in Nxrange:
    for Ny in Nyrange:
        globalstrain = laminate1.globalstrains([Nx, Ny, 0, 0, 0, 0])
        loadcombinations.update({(Nx, Ny): False})

        for ply in laminaarray:
            plystrainsG = ply.plystrains(globalstrain, 1)
            plystressesG = ply.calculatestresses(plystrainsG)
            plystressesL = ply.principalstresses2(plystressesG)
            sigma11, sigma22, sigma12 = plystressesL
            print(plystressesL)

            # dFF = ply.puckFF(Xt, Xc, E1, 200E3, v12, 0.1, sigma11, sigma22, sigma12)
            # dIFFA = ply.puckIFFA(sigma12, sigma22, S12, Yt)
            # dIFFB = ply.puckIFFB(sigma12, sigma22, S12, Yc)
            # dIFFC = ply.puckIFFC(sigma12, sigma22, Yc, S12)
            # dlist = [dFF, dIFFA, dIFFB, dIFFC]

            d1 = ply.hashinFT(Xt, S12, sigma11, sigma12)
            d2 = ply.hashinFC(Xc, sigma11)
            d3 = ply.hashinMT(Yt, S12, sigma22, sigma12)
            d4 = ply.hashinMC(Yc, S12, sigma22, sigma12)
            dlist = [d1, d2, d3, d4]
            d = max(dlist)
            print(d, Nx, Ny)
            if d >= 1:
                loadcombinations.update({(Nx, Ny): True})


for key in loadcombinations.keys():
    if loadcombinations[key] is True:
        envelopepoints.append(np.array(key))

envelopepoints = np.array(envelopepoints)

plt.scatter(envelopepoints[:, 0], envelopepoints[:, 1])
plt.ylabel("Ny (N/mm)")
plt.xlabel("Nx (N/mm)")
plt.xlim(-1000, 1000)
plt.ylim(-1000, 1000)
plt.grid()
plt.show()

print(envelopepoints)



























