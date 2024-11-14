""" MODULE TO GET RESONANCE LOCATIONS AND FREQUENCIES """

import numpy as np

def r_LR(m, l, ab=1, outer=True):
    """ LINDBLAD RESONANCES """
    if outer:
        eps = 1
    else:
        eps = -1
    return(((m + eps) / (l))**(2/3.) * ab)

def r_CR(m, l, ab=1):
    """ CO-ROTATION RESONANCES """
    return((m/l)**(2./3.)*ab)

def omega(r, mb=1, G=1):
    return(np.sqrt(G*mb/(r**3)))

def r_from_omega(om, mb=1, G=1):
    return((G*mb/(om**2))**(1/3))

