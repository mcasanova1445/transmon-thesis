
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *
import tgates

def R2(psi0, c1, c2):                 # 1: black; 0: white; -1: none
    if c1 == -1 and c2 == -1:
        res = tgates.CNOT(psi0, 3, 2)
        res = tgates.X(res.states[-1], 3)

    elif c1 == -1:
        psi1 = psi0

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

        res = tgates.CCNOT(ps1, 1, 3, 2)
        res = tgates.CNOT(res.states[-1], 1, 3)
        psi1 = tgates.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    elif c2 == -1:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        res = tgates.CCNOT(ps1, 0, 3, 2)
        res = tgates.CNOT(res.states[-1], 0, 3)
        psi1 = tgates.states[-1]

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

    else:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

        res = tgates.CCCNOT(ps1, 0, 1, 3, 2)
        res = tgates.CCNOT(res.states[-1], 0, 1, 3)
        psi1 = tgates.states[-1]

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    return res

def L2(psi0, c1, c2):                 # 1: black; 0: white; -1: none
    if c1 == -1 and c2 == -1:
        res = tgates.X(psi0, 3)
        res = tgates.CNOT(res.states[-1], 3, 2)

    elif c1 == -1:
        psi1 = psi0

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

        res = tgates.CNOT(ps1, 1, 3)
        res = tgates.CCNOT(res.states[-1], 1, 3, 2)
        psi1 = tgates.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    elif c2 == -1:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        res = tgates.CNOT(ps1, 0, 3)
        res = tgates.CCNOT(res.states[-1], 0, 3, 2)
        psi1 = tgates.states[-1]

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

    else:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

        res = tgates.CCNOT(ps1, 0, 1, 3)
        res = tgates.CCCNOT(res.states[-1], 0, 1, 3, 2)
        psi1 = tgates.states[-1]

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    return res

def T1(psi0):
    res = L2(psi0, 0, 1)
    return L2(res.states[-1], 1, 0)

def Kb2(psi0):                          #CCH
    res = CCRy(psi0, 0, 1, 2, np.pi/2)
    res = CCNOT(res.states[-1], 0, 1, 2)
    res = CCRy(res.states[-1], 0, 1, 3, np.pi/2)
    return CCNOT(res.states[-1], 0, 1, 3)

def Kb1(psi0):
    thetay00 = -1.900220607202805
    thetay01 = -2.212560440981948
    thetay10 = -0.75 * np.pi

    res = X(psi0, 3)
    res = CRy(res.states[-1], 3, 2, thetay01)
    res = X(res.states[-1], 3)
    res = X(res.states[-1], 2)
    res = CRy(res.states[-1], 2, 3, thetay00)
    res = X(res.states[-1], 2)
    res = CRy(res.states[-1], 2, 3, thetay10)

    res = CCCRy(res.states[-1], 0, 1, 2, 3, -thetay10)
    res = CCNOT(res.states[-1], 0, 1, 2)
    res = CCCRy(res.states[-1], 0, 1, 2, 3, -thetay00)
    res = CCNOT(res.states[-1], 0, 1, 2)
    res = CCNOT(res.states[-1], 0, 1, 3)
    res = CCCRy(res.states[-1], 0, 1, 3, 2, -thetay01)
    return CCNOT(res.states[-1], 0, 1, 3)

def Dd(psi0):
    return 
