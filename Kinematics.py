import ROOT
from ROOT import TMath
import numpy as np

def eta(p):
    pt = np.sqrt(p[1]**2+p[2]**2)
    pz = p[3]
    return np.log(abs(np.tan((pt-pz)/(2*(pt+pz)))))

def photonpt(d):
    return np.sqrt(d[1]**2 + d[2]**2)

def invariantmass(d):
    return np.sqrt((d[0]+d[4])**2-(d[1]+d[5])**2-(d[2]+d[6])**2-(d[3]+d[7])**2)

def TInvariant(d):
    return np.sqrt(abs((d[0]-d[4])**2-(d[1]-d[5])**2-(d[2]-d[6])**2-(d[3]-d[7])**2))
