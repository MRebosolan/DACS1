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

Exarray = np.zeros(10)    #EACH ROW IS A COLUMN FOR A SPECIFIC THETA, n ROWS
Eyarray = np.zeros(10)
vxyarray = np.zeros(10)
Gxy = np.zeros(10)
Efxarray = np.zeros(10)
Efyarray = np.zeros(10)
vfxyarray = np.zeros(10)
Gfxy = np.zeros(10)


for n in range(1, 11):

    theta_array = np.linspace(-90, 90, 181)



    for theta in theta_array:
        angles = []
        for i in range(n):

            angles.append(30)
            angles.append(theta)
            angles.append(-theta)
            angles.append(60)
            angles.append(60)


        Exlist = []  #ARRAYS FOR ONE VALUE OF N
        Eylist = []
        vxylist = []
        Gxylist = []
        Efxlist = []
        Efylist = []
        vfxylist = []
        Gfxylist = []
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
        Efylsit = np.array(Efylist)
        vfxylist = np.array(vfxylist)
        Gfxylist = np.array(Gfxylist)
        print(n)




