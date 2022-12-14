from Functions import *
from math import *
import numpy as np
import matplotlib.pyplot as plt


E1 = 140.03E3
E2 = 8.036E3
v12 = 0.355
G12 = 4.363E3
S12 = 98.92
Xt = 1926.1
Xc = 1480
Yt = 107.06
Yc = 220


laminaarray = []
angles = []


angles.append(30)
angles.append(45)
angles.append(-45)
angles.append(60)
angles.append(60)


t = 2/len(angles)
zlocations = np.arange(-t*len(angles)/2, t*(len(angles)+1)/2, t)



for i in range(len(angles)):
    ply = lamina(E1, E2, v12, G12, radians(angles[i]), zlocations[i], zlocations[i+1], Xt, Xc, Yt, Yc, S12)
    laminaarray.append(ply)

laminate1 = laminate(laminaarray)
ABD = laminate1.createABD()
ABDI = np.linalg.inv(ABD)

loads = [20, 5, 0, 0, 9, 0]
globalstrain = laminate1.globalstrains(loads)


strainarray = []
sigmaarray = []
globalstrainarray = []

for ply in laminaarray:
    plystrainG = ply.plystrains(globalstrain, 10)

    for point in range(10):
        plystrainGatpoint = plystrainG[:, point]
        plystressGatpoint = ply.calculatestresses(plystrainGatpoint)
        plystrainLatpoint = ply.principalstrains(plystrainGatpoint)
        plystressLatpoint = ply.principalstresses2(plystressGatpoint)
        globalstrainarray.append(plystrainGatpoint)
        strainarray.append(plystrainLatpoint)
        sigmaarray.append(plystressLatpoint)


strainarray = np.array(strainarray)
sigmaarray = np.array(sigmaarray)
globalstrainarray = np.array(globalstrainarray)


plt.plot(strainarray[:, 0], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "r", label="strain11")
plt.plot(strainarray[:, 1], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "g", label="strain22")
plt.plot(strainarray[:, 2], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "b", label="strain12")

# plt.plot(sigmaarray[:, 0], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "r", label="sigma11")
# plt.plot(sigmaarray[:, 1], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "g", label="sigma22")
# plt.plot(sigmaarray[:, 2], np.linspace(zlocations[0], zlocations[-1], len(sigmaarray[:, 0])), "b", label="sigma12")


plt.grid()
plt.xlabel("strain (-)")
plt.ylabel("z (mm)")
plt.legend()
plt.show()
