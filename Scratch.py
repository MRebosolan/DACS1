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
v12 = 0.35
G12 = 4.685E3
S12 = 79
S23 = 20
Xt = 1950
Xc = 1480
Yt = 107
Yc = 220

t = 0.125
puck = True


zlocations = np.arange(-t * len(angles) / 2, t * (len(angles) + 1) / 2, t)

laminaarray = []

for n in range(len(angles)):
    ply = lamina(E1, E2, v12, G12, radians(angles[n]), zlocations[n], zlocations[n + 1], Xt, Xc, Yt, Yc, S12)
    laminaarray.append(ply)

laminate1 = laminate(laminaarray)

Nxrange = np.arange(-2000, 2000, 25)
Nyrange = np.arange(0, -2000, -25)
Nyrange2 = np.arange(0, 2000, 25)

envelopepoints = []
envelopepoints2 = []
loadcombinations = {}  # KEYS = ENVELOPE COORDINATES, VALUES = FAILED TRUE/FALSE
loadcombinations2 = {}
for Nx in Nxrange:
    for Ny in Nyrange2:
        loadcombinations2.update({(Nx, Ny): False})


# PLYSTRAINS --> CALCULATESTRESSES --> PRINCIPALSTRESSES2

# for Nx in Nxrange:
#     for ply in laminaarray:
#         ply.failed = False
#     for Ny in Nyrange:
#         globalstrain = laminate1.globalstrains([Nx, Ny, 0, 0, 0, 0])
#         loadcombinations.update({(Nx, Ny): False})
#
#         for ply in laminaarray:
#             plystrainsG = ply.plystrains(globalstrain, 1)
#             plystressesG = ply.calculatestresses(plystrainsG)
#             plystressesL = ply.principalstresses2(plystressesG)
#             sigma11, sigma22, sigma12 = plystressesL
#
#             if puck is True:
#                 dFF = ply.puckFF(Xt, Xc, E1, 200E3, v12, 0.1, sigma11, sigma22, sigma12)
#                 dIFFA = ply.puckIFFA(sigma12, sigma22, S12, Yt)
#                 dIFFB = ply.puckIFFB(sigma12, sigma22, S12, Yc)
#                 dIFFC = ply.puckIFFC(sigma12, sigma22, Yc, S12)
#                 dlist = [dFF, dIFFA, dIFFB]#, dIFFC]
#             else:
#                 d1 = ply.hashinFT(Xt, S12, sigma11, sigma12)
#                 d2 = ply.hashinFC(Xc, sigma11)
#                 d3 = ply.hashinMT(Yt, S12, sigma22, sigma12)
#                 d4 = ply.hashinMC(Yc, S12, S23, sigma22, sigma12)
#                 dlist = [d1, d2, d3, d4]
#             d = max(dlist)
#
#             if d >= 1:
#                 loadcombinations.update({(Nx, Ny): True})
#                 ply.failed = True
#         if any(ply.failed is True for ply in laminaarray):
#             break
#
# for key in loadcombinations.keys():
#     if loadcombinations[key] is True:
#         envelopepoints.append(np.array(key))


for Nx in Nxrange:
    for ply in laminaarray:
        ply.failed = False
        ply.dtype = 4
    for Ny in Nyrange2:
        for ply in laminaarray:
            if ply.failed is True:
                ply.E1 = 0.001
                ply.E2 = 0.001
                ply.Xt = 0.001
                ply.Xc = 0.001
                ply.Yt = 0.001
                ply.Yt = 0.001
                ply.S12 = 0.001
                ply.G12 = 0.001

        globalstrain = laminate1.globalstrains([Nx, Ny, 0, 0, 0, 0])

        for ply in laminaarray:
            plystrainsG = ply.plystrains(globalstrain, 1)
            plystressesG = ply.calculatestresses(plystrainsG)
            plystressesL = ply.principalstresses2(plystressesG)
            sigma11, sigma22, sigma12 = plystressesL

            if puck is True:
                dFF = ply.puckFF(Xt, Xc, E1, 200E3, v12, 0.1, sigma11, sigma22, sigma12)
                dIFFA = ply.puckIFFA(sigma12, sigma22, S12, Yt)
                dIFFB = ply.puckIFFB(sigma12, sigma22, S12, Yc)
                dIFFC = ply.puckIFFC(sigma12, sigma22, Yc, S12)
                dlist = [dFF, dIFFA, dIFFB]#, dIFFC]
            else:
                d1 = ply.hashinFT(Xt, S12, sigma11, sigma12)
                d2 = ply.hashinFC(Xc, sigma11)
                d3 = ply.hashinMT(Yt, S12, sigma22, sigma12)
                d4 = ply.hashinMC(Yc, S12, S23, sigma22, sigma12)
                dlist = [d1, d2, d3, d4]

            d = max(dlist)
            dtype = dlist.index(d)


            if d >= 1 and ply.failed is False:
                loadcombinations2.update({(Nx, Ny): True})
                ply.failed = True
                ply.dtype = dtype
            print(dlist, Nx, Ny, ply.dtype, ply.failed)

        if any(ply.failed is True for ply in laminaarray):
            break

for key in loadcombinations2.keys():
    if loadcombinations2[key] is True:
        envelopepoints2.append(np.array(key))

envelopepoints = np.array(envelopepoints)
envelopepoints2 = np.array(envelopepoints2)

#plt.plot(envelopepoints[:, 0], envelopepoints[:, 1], "b")
plt.plot(envelopepoints2[:, 0], envelopepoints2[:, 1], "b")
plt.ylabel("Ny (N/mm)")
plt.xlabel("Nx (N/mm)")
plt.xlim(-2000, 2000)
plt.ylim(-2000, 2000)
plt.grid()
plt.show()



