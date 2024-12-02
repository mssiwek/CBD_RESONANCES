""" MODULE TO GET RESONANCE LOCATIONS AND FREQUENCIES """

import numpy as np

def r_ILR(l,m, ab=1):
    """ INNER LINDBLAD RESONANCES """
    eps = -1
    if l == 0 or m <= 1:
        return(np.nan)
    else:
        return(((m + eps) / (l))**(2/3.) * ab)

def r_OLR(l,m, ab=1):
    """ OUTER LINDBLAD RESONANCES """
    eps = +1
    if l == 0 or m == -1:
        return(np.nan)
    else:
        return(((m + eps) / (l))**(2/3.) * ab)


def r_CR(l,m, ab=1):
    """ CO-ROTATION RESONANCES """
    if l == 0 or m == 0:
        return(np.nan)
    else:
        return((m/l)**(2./3.)*ab)

def omega(r, mb=1, G=1):
    return(np.sqrt(G*mb/(r**3)))

def r_from_omega(om, mb=1, G=1):
    return((G*mb/(om**2))**(1/3))

