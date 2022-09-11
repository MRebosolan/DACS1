
def hashinMT(Yt, S12, sigma2, tau12):
        if sigma2 > 0:
            d = (sigma2**2)/(Yt**2) + (tau12**2)/(S12**2)
        else:
            d = 0
        return d

def hashinMC(Yc, S12, S23, sigma2, tau12):
    if sigma2 < 0:
        d = (sigma2/Yc)*((Yc/(2*S23))**2 - 1) + (sigma2/(2*S23))**2 + (tau12/S12)**2
    else:
        d = 0
    return d

print(hashinMT(107, 79, 50, 25))
print(hashinMC(220, 79, 300, -50, 25))