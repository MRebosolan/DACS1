import numpy as np
from math import *
import matplotlib.pyplot as plt
from Functions import *
import pandas as pd

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
E2 = 8.036E3
v12 = 0.355
G12 = 4.363E3
S12 = 98.92
S23 = 20
Xt = 1926.1
Xc = 1480
Yt = 107.06
Yc = 220


t = 2/len(angles)
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
LPFcoordinates = {}
LPFenvelope = []
LPFcoordinates2 = {}
LPFenvelope2 = []
LPFstrains = []

for Nx in Nxrange:
    for Ny2 in Nyrange2:
        loadcombinations2.update({(Nx, Ny2): False})
        LPFcoordinates2.update({(Nx, Ny2): False})
    for Ny in Nyrange:
        loadcombinations.update({(Nx, Ny): False})
        LPFcoordinates.update({(Nx, Ny): False})


#PLYSTRAINS --> CALCULATESTRESSES --> PRINCIPALSTRESSES2

for Nx in Nxrange:
    FPF = False
    for ply in laminaarray:
        ply.failed = False
        ply.FF = False
        ply.MF = False
        ply.dtype = 4
    for Ny in Nyrange:
        for ply in laminaarray:
            if ply.failed is True:
                if ply.FF is True:
                    ply.E1 = 0.001
                    ply.E2 = 0.001
                    ply.Xt = 0.001
                    ply.Xc = 0.001
                    ply.Yt = 0.001
                    ply.Yc = 0.001
                    ply.S12 = 0.001
                    ply.G12 = 0.001
                if ply.MF is True:
                    ply.E2 = 0.85 * E2
                    ply.Yt = 0.85 * Yt
                    ply.Yc = 0.85 * Yc
                    ply.S12 = 0.85 * S12
                    ply.G12 = 0.85 * G12
            else:
                ply.E1 = E1
                ply.E2 = E2
                ply.Xt = Xt
                ply.Xc = Xc
                ply.Yt = Yt
                ply.Yc = Yc
                ply.S12 = S12
                ply.G12 = G12

        globalstrain = laminate1.globalstrains([Nx, Ny, 0, 0, 0, 0])

        for ply in laminaarray:
            plystrainsG = ply.plystrains(globalstrain, 1)
            plystressesG = ply.calculatestresses(plystrainsG)
            plystressesL = ply.principalstresses2(plystressesG)
            sigma11, sigma22, sigma12 = plystressesL

            if puck is True:
                dFF = ply.puckFF(Xt, Xc, E1, 225E3, v12, 0.2, sigma11, sigma22, sigma12)
                dIFFA = ply.puckIFFA(sigma12, sigma22, S12, Yt)
                dIFFB = ply.puckIFFB(sigma12, sigma22, S12, Yc)
                dIFFC = ply.puckIFFC(sigma12, sigma22, Yc, S12)
                dlist = [dFF, dIFFA, dIFFB, dIFFC]
            else:
                d1 = ply.hashinFT(sigma11, sigma12)
                d2 = ply.hashinFC(sigma11)
                d3 = ply.hashinMT(sigma22, sigma12)
                d4 = ply.hashinMC(sigma22, sigma12)
                dlist = [d1, d2, d3, d4]

            d = max(dlist)
            dtype = dlist.index(d)


            if d >= 1 and ply.failed is False:
                ply.failed = True
                if dlist[0] >= 1:
                    ply.FF = True
                else:
                    ply.MF = True
                if FPF is False:
                    loadcombinations.update({(Nx, Ny): True})
                FPF = True

            #print(dlist, Nx, Ny, degrees(ply.theta), ply.failed, FPF)

        if all(ply.failed is True for ply in laminaarray):
            LPFcoordinates.update({(Nx, Ny): True})
            LPFstrains.append([Nx, Ny, *globalstrain[:3]])
            break

for key in loadcombinations.keys():
    if loadcombinations[key] is True:
        envelopepoints.append(np.array(key))

for key in LPFcoordinates.keys():
    if LPFcoordinates[key] is True:
        LPFenvelope.append(np.array(key))



for Nx in Nxrange:
    FPF2 = False
    for ply in laminaarray:
        ply.failed = False
        ply.FF = False
        ply.MF = False
        ply.dtype = 4
    for Ny in Nyrange2:
        for ply in laminaarray:
            if ply.failed is True:
                if ply.FF is True:
                    ply.E1 = 0.001
                    ply.E2 = 0.001
                    ply.Xt = 0.001
                    ply.Xc = 0.001
                    ply.Yt = 0.001
                    ply.Yc = 0.001
                    ply.S12 = 0.001
                    ply.G12 = 0.001
                if ply.MF is True:
                    ply.E2 = 0.85 * E2
                    ply.Yt = 0.85 * Yt
                    ply.Yc = 0.85 * Yc
                    ply.S12 = 0.85 * S12
                    ply.G12 = 0.85 * G12
            else:
                ply.E1 = E1
                ply.E2 = E2
                ply.Xt = Xt
                ply.Xc = Xc
                ply.Yt = Yt
                ply.Yc = Yc
                ply.S12 = S12
                ply.G12 = G12

        globalstrain = laminate1.globalstrains([Nx, Ny, 0, 0, 0, 0])

        for ply in laminaarray:
            plystrainsG = ply.plystrains(globalstrain, 1)
            plystressesG = ply.calculatestresses(plystrainsG)
            plystressesL = ply.principalstresses2(plystressesG)
            sigma11, sigma22, sigma12 = plystressesL

            if puck is True:
                dFF = ply.puckFF(Xt, Xc, E1, 225E3, v12, 0.2, sigma11, sigma22, sigma12)
                dIFFA = ply.puckIFFA(sigma12, sigma22, S12, Yt)
                dIFFB = ply.puckIFFB(sigma12, sigma22, S12, Yc)
                dIFFC = ply.puckIFFC(sigma12, sigma22, Yc, S12)
                dlist = [dFF, dIFFA, dIFFB, dIFFC]
            else:
                d1 = ply.hashinFT(sigma11, sigma12)
                d2 = ply.hashinFC(sigma11)
                d3 = ply.hashinMT(sigma22, sigma12)
                d4 = ply.hashinMC(sigma22, sigma12)
                dlist = [d1, d2, d3, d4]

            d = max(dlist)
            dtype = dlist.index(d)


            if d >= 1 and ply.failed is False:
                ply.failed = True
                if dlist[0] >= 1:
                    ply.FF = True
                else:
                    ply.MF = True
                if FPF2 is False:
                    loadcombinations2.update({(Nx, Ny): True})
                FPF2 = True

            print(dlist, Nx, Ny, dtype, degrees(ply.theta), ply.failed, FPF2)

        if all(ply.failed is True for ply in laminaarray):
            LPFcoordinates2.update({(Nx, Ny): True})
            LPFstrains.append([Nx, Ny, *globalstrain[:3]])
            break

for key in loadcombinations2.keys():
    if loadcombinations2[key] is True:
        envelopepoints2.append(np.array(key))

for key in LPFcoordinates2.keys():
    if LPFcoordinates2[key] is True:
        LPFenvelope2.append(np.array(key))

envelopepoints = np.array(envelopepoints)
envelopepoints2 = np.array(envelopepoints2)
LPFenvelope = np.array(LPFenvelope)
LPFenvelope2 = np.array(LPFenvelope2)
LPFstrains = np.array(LPFstrains)


df = pd.DataFrame()
df["Nx"] = LPFstrains[:, 0]
df["Ny"] = LPFstrains[:, 1]
df["Exx"] = LPFstrains[:, 2]
df["Eyy"] = LPFstrains[:, 3]
df["Exy"] = LPFstrains[:, 4]

filename = "GlobalstrainsNxNywithout.xlsx"
df.to_excel(filename)

plt.plot(envelopepoints[:, 0], envelopepoints[:, 1], "b")
plt.plot(LPFenvelope[:, 0], LPFenvelope[:, 1], "r")
plt.plot(envelopepoints2[:, 0], envelopepoints2[:, 1], "b")
plt.plot(LPFenvelope2[:, 0], LPFenvelope2[:, 1], "r")
plt.ylabel("Ny (N/mm)")
plt.xlabel("Nx (N/mm)")
plt.xlim(-2000, 2000)
plt.ylim(-2000, 2000)
plt.grid()
plt.show()