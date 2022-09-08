from Functions import *
from math import *
import numpy as np
import matplotlib.pyplot as plt


E1 = 140.03E3
E2 = 7.72E3
v12 = 0.3
G12 = 4.685E3
S12 = 79
Xt = 1950
Xc = 1200
Yt = 48
Yc = 35

Exarray = np.zeros((10, 181))
Eyarray = np.zeros((10, 181))
vxyarray = np.zeros((10, 181))
Gxyarray = np.zeros((10, 181))
Efxarray = np.zeros((10, 181))
Efyarray = np.zeros((10, 181))
vfxyarray = np.zeros((10, 181))
Gfxyarray = np.zeros((10, 181))

theta_array = np.linspace(-90, 90, 181)


for n in range(1, 11):



    Exlist = []  # ARRAYS FOR ONE VALUE OF N
    Eylist = []
    vxylist = []
    Gxylist = []
    Efxlist = []
    Efylist = []
    vfxylist = []
    Gfxylist = []


    for theta in theta_array:
        angles = []
        for i in range(n):

            angles.append(30)
            angles.append(theta)
            angles.append(-theta)
            angles.append(60)
            angles.append(60)



        t = 0.125
        zlocations = np.arange(-t*len(angles)/2, t*(len(angles)+1)/2, t)

        laminaarray = []

        for j in range(len(angles)):
            ply = lamina(E1, E2, v12, G12, radians(angles[j]), zlocations[j], zlocations[j + 1], Xt, Xc, Yt, Yc, S12)
            laminaarray.append(ply)

        laminate1 = laminate(laminaarray)
        Ex, Ey, vxy, Gxy, Efx, Efy, vfxy, Gfxy = laminate1.equivalentelastic()
        Exlist.append(Ex)
        Eylist.append(Ey)
        vxylist.append(vxy)
        Gxylist.append(Gxy)
        Efxlist.append(Efx)
        Efylist.append(Efy)
        vfxylist.append(vfxy)
        Gfxylist.append(Gfxy)



    Exlist = np.array(Exlist)
    Eylist = np.array(Eylist)
    vxylist = np.array(vxylist)
    Gxylist = np.array(Gxylist)
    Efxlist = np.array(Efxlist)
    Efylist = np.array(Efylist)
    vfxylist = np.array(vfxylist)
    Gfxylist = np.array(Gfxylist)

    Exarray[n-1] = Exlist
    Eyarray[n-1] = Eylist
    vxyarray[n-1] = vxylist
    Gxyarray[n-1] = Gxylist
    Efxarray[n-1] = Efxlist
    Efyarray[n-1] = Efylist
    vfxyarray[n-1] = vfxylist
    Gfxyarray[n-1] = Gfxylist

for q, row in enumerate(vfxyarray):
    plt.plot(theta_array, row, label=f"n={q+1}")

plt.xlabel("theta (deg)")
plt.ylabel("vfxy (-)")
plt.xlim(-89, 90)
plt.grid()
plt.legend()
plt.show()











