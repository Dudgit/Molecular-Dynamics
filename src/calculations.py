import numpy as np
from spce_vars import *

def LJ(eps, sig, r):
    return (4*eps*((sig/r)**12-(sig/r)**6))

def Coul(qi, qj, r):
    #eps_0 = 8.8541878128*10**(-12) farad/meter = 8.8541878128*10**(-2) farad/angstrom
    return (qi*qj)/(4*np.pi*(8.8541878128*10**(-2))*r)

def LJC(eps, sig, qi, qj, r):
    return LJ(eps, sig, r) + Coul(qi, qj, r)

def LBcomb(epsi, epsj, sigi, sigj):
    epsij = (epsi*epsj)**0.5
    sigij = (sigi+sigj)/2
    return (epsij, sigij)

def dFH2(dFH1):
    return (dFH1**2 + 2.6666*dFH1 + 2.6666)**0.5