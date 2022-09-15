import numpy as np
from math import *



class lamina:
    def __init__(self, E1, E2, v12, G12, theta, z1, z2, Xt, Xc, Yt, Yc, S12, failed=False, dtype=4, FF=False, MF=False):
        self.failed = failed
        self.dtype = dtype
        self.E1 = E1
        self.E2 = E2
        self.v12 = v12
        self.v21 = v12 * (E2/E1)
        self.G12 = G12
        self.theta = theta
        self.z1 = z1
        self.z2 = z2
        self.Xt = Xt
        self.Xc = Xc
        self.Yt = Yt
        self.Yc = Yc
        self.S12 = S12
        self.FF = FF
        self.MF = MF


    def Qmatrix(self):
        Q = 1 - self.v12 * self.v21
        Q11 = self.E1 / Q
        Q22 = self.E2 / Q
        Q12 = self.v12 * self.E2 / Q
        Q66 = self.G12

        m = cos(self.theta)
        n = sin(self.theta)

        Qxx = Q11 * m ** 4 + 2 * (Q12 + 2 * Q66) * m * m * n * n + Q22 * n ** 4
        Qxy = (Q11 + Q22 - 4 * Q66) * m * m * n * n + Q12*(m ** 4 + n ** 4)
        Qyy = Q11 * n ** 4 + 2 * (Q12 + 2 * Q66) * m * m * n * n + Q22 * m ** 4
        Qxs = (Q11 - Q12 - 2 * Q66) * n * m ** 3 + (Q12 - Q22 + 2 * Q66) * n ** 3 * m
        Qys = (Q11 - Q12 - 2 * Q66) * m * n ** 3 + (Q12 - Q22 + 2 * Q66) * m ** 3 * n
        Qss = (Q11 + Q22 - 2 * Q12 - 2 * Q66) * n * n * m * m + Q66 * (n ** 4 + m ** 4)

        Qmatrix = np.array([[Qxx, Qxy, Qxs], [Qxy, Qyy, Qys], [Qxs, Qys, Qss]]) #lamina Q matrix in global coordinates

        return Qmatrix

    def plystrains(self, globalstrains, n):            #CONVERTS LAMINATE STRAINS (ABD OUTPUT) TO LAMINA STRAINS
        exx_bottom = globalstrains[0] + globalstrains[3] * self.z1
        exx_top = globalstrains[0] + globalstrains[3] * self.z2
        eyy_bottom = globalstrains[1] + globalstrains[4] * self.z1
        eyy_top = globalstrains[1] + globalstrains[4] * self.z2
        es_bottom = globalstrains[2] + globalstrains[5] * self.z1
        es_top = globalstrains[2] + globalstrains[5] * self.z2
        if n == 1:
            exx = (exx_top + exx_bottom)/2
            eyy = (eyy_top + eyy_bottom)/2
            es = (es_top + es_bottom)/2
        else:
            exx = np.linspace(exx_bottom, exx_top, n)
            eyy = np.linspace(eyy_bottom, eyy_top, n)
            es = np.linspace(es_bottom, es_top, n)
        return np.array([exx, eyy, es])

    def principalstrains(self, strains):   #ROTATES PLY STRAINS FROM GLOBAL CS TO PRINCIPAL CS

        m = cos(self.theta)
        n = sin(self.theta)

        rm = np.array([[m**2, n**2, m*n], [n**2, m**2, -m*n], [-2*m*n, 2*m*n, m**2 - n**2]])
        principalstrains = rm @ strains
        return principalstrains


    def calculatestresses(self, globalstrains):     #stresses in global coordinates from strains in global coordinates
        Qmat = self.Qmatrix()
        stresses = Qmat @ globalstrains
        return stresses

    def principalstresses(self, principalstrains):


        Q = 1 - self.v12 * self.v21
        Q11 = self.E1 / Q
        Q22 = self.E2 / Q
        Q12 = self.v12 * self.E2 / Q
        Q66 = self.G12
        Qmat = np.array([[Q11, Q12, 0], [Q12, Q22, 0], [0, 0, Q66]])
        principalstresses = Qmat @ principalstrains
        return principalstresses

    def principalstresses2(self, globalstresses): #ROTATES PLY STRESSES FROM GLOBAL CS TO PRINCIPAL CS

        m = cos(self.theta)
        n = sin(self.theta)

        rm = np.array([[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, m**2-n**2]])
        principalstresses = rm @globalstresses
        return principalstresses


    def hashinFT(self, sigma1, tau12):
        Xt = self.Xt
        S12 = self.S12
        if sigma1 >= 0:
            d = (sigma1/Xt)**2 + (tau12/S12)**2
        else:
            d = 0
        return d

    def hashinFC(self, sigma1):
        Xc = self.Xc
        if sigma1 < 0:
            d = -sigma1/Xc
        else:
            d = 0
        return d

    def hashinMT(self, sigma2, tau12):
        Yt = self.Yt
        S12 = self.S12
        if sigma2 > 0:
            d = (sigma2**2)/(Yt**2) + (tau12**2)/(S12**2)
        else:
            d = 0
        return d

    def hashinMC(self, sigma2, tau12, S23 = 20):
        Yc = self.Yc
        S12 = self.S12
        if sigma2 < 0:
            d = (sigma2/Yc)*((Yc/(2*S23))**2 - 1) + (sigma2/(2*S23))**2 + (tau12/S12)**2
        else:
            d = 0
        return d



    def puckFF(self, Xt, Xc, E1, E1f, v12, v12f, sigma1, sigma2, mf=1.3):  #mf=1.3 for GFRP and 1.1 for CFRP
        if sigma1 > 0:
            d = (1/Xt)*(sigma1 - (v12 - v12f*mf*(E1/E1f))*sigma2)
        else:
            d = (-1 / Xc) * (sigma1 - (v12 - v12f * mf * (E1 / E1f)) * sigma2)
        return d

    def puckIFFA(self, sigma21, sigma22, S21, Yt, p12=0.3):
        if sigma22 > 0:
            d = sqrt((sigma21/S21)**2 + (1-p12*(Yt/S21))**2 * (sigma22/Yt)**2) + p12*(sigma22/S21)   #p12=0.3
        else:
            d = 0
        return d

    def puckIFFB(self,  sigma12, sigma22, S12, Yc, p12=0.2, p23=0.25):

        sigma23a = (S12/(2*p12))*(sqrt(1+2*p12*(Yc/S12)) - 1)
        sigma12c = S12*sqrt(1+2*p23)
        if 0 <= abs(sigma22 / sigma12) <= sigma23a/abs(sigma12c) and sigma22 < 0:
            #d = 1/S12 * ((sqrt(sigma12**2 + (p12*sigma22)**2)) + p12*sigma22)  #p12 = 0.2
            d = sqrt((sigma12/S12)**2 + ((p12/S12)*sigma22)**2) + (p12/S12)*sigma22
        else:
            d = 0
        return d

    def puckIFFC(self, sigma12, sigma22, Yc, S12, p12=0.2, p23=0.25):   #p3 = 0.25
        sigma23a = (S12 / (2 * p12)) * (sqrt(1 + 2 * p12 * (Yc / S12)) - 1)
        sigma12c = S12 * sqrt(1 + 2 * p23)

        if 0 <= abs(sigma12/sigma22) <= abs(sigma12c)/sigma23a and sigma22 < 0:
            d = (((sigma12/(2*(1+p23)*S12)))**2 + (sigma22/Yc)**2)*(-Yc/sigma22)
        else:
            d = 0
        return d



class laminate:

    def __init__(self, laminaarray):
        self.laminaarray = laminaarray
        self.nplies = len(laminaarray)
        self.h = laminaarray[-1].z2 - laminaarray[0].z1



    def createABD(self):
        A = np.zeros((3, 3))
        B = np.zeros((3, 3))
        D = np.zeros((3, 3))
        for k, ply in enumerate(self.laminaarray):
            A += ply.Qmatrix() * (ply.z2 - ply.z1)
            B += 0.5*(ply.Qmatrix() * (ply.z2**2 - ply.z1**2))
            D += 1/3 * (ply.Qmatrix()) * (ply.z2**3 - ply.z1**3)


        AB = np.hstack((A, B))
        BD = np.hstack((B, D))
        ABD = np.vstack((AB, BD))
        return ABD

    def globalstrains(self, loads):
        ABD = self.createABD()
        loads = np.array(loads)
        strains = np.linalg.inv(ABD) @ loads
        return strains


    def equivalentelastic(self):
        ABD = self.createABD()
        ABDI = np.linalg.inv(ABD)
        h = self.h

        Ex = 1/(ABDI[0, 0] * h)
        Ey = 1/(ABDI[1, 1] * h)
        vxy = -(ABDI[0, 1])/(ABDI[0, 0])
        Gxy = 1/(ABDI[2, 2] * h)
        Efx = 12/(ABDI[3, 3] * h**3)
        Efy = 12/(ABDI[4, 4] * h**3)
        vfxy = - ABDI[3, 4]/ABDI[3, 3]
        Gfxy = 12/(ABDI[5, 5] * h**3)

        return [Ex, Ey, vxy, Gxy, Efx, Efy, vfxy, Gfxy]


















