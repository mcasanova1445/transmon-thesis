# PageRank del grafo estrella

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *
import tgates
import time

def R2(psi0, c1, c2):                 # 1: black; 0: white; -1: none
    if c1 == -1 and c2 == -1:
        res = tgates.CNOT(psi0, 3, 2)
        res = tgates.X(res.states[-1], 3)

    elif c1 == -1:
        psi1 = psi0

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

        res = tgates.CCNOT(psi1, 1, 3, 2)
        res = tgates.CNOT(res.states[-1], 1, 3)
        psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    elif c2 == -1:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        res = tgates.CCNOT(psi1, 0, 3, 2)
        res = tgates.CNOT(res.states[-1], 0, 3)
        psi1 = res.states[-1]

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

        res = tgates.CCCNOT(psi1, 0, 1, 3, 2)
        res = tgates.CCNOT(res.states[-1], 0, 1, 3)
        psi1 = res.states[-1]

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

        res = tgates.CNOT(psi1, 1, 3)
        res = tgates.CCNOT(res.states[-1], 1, 3, 2)
        psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    elif c2 == -1:
        psi1 = psi0

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        res = tgates.CNOT(psi1, 0, 3)
        res = tgates.CCNOT(res.states[-1], 0, 3, 2)
        psi1 = res.states[-1]

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

        res = tgates.CCNOT(psi1, 0, 1, 3)
        res = tgates.CCCNOT(res.states[-1], 0, 1, 3, 2)
        psi1 = res.states[-1]

        if c1 == 0:
            res = tgates.X(psi1, 0)
            psi1 = res.states[-1]

        if c2 == 0:
            res = tgates.X(psi1, 1)
            psi1 = res.states[-1]

    return res

def Kb1(psi0):
    # theta target trigger_state
    # \[Theta]y00 -> 1.85806, \[Theta]y10 -> \ 2.48274, \[Theta]y11 -> 1.5708
    thetay00 = 1.85806
    thetay10 = 2.48274
    thetay11 = 1.5708

    res = tgates.X(psi0, 0)
    res = tgates.X(res.states[-1], 1)

    res = tgates.X(res.states[-1], 3)
    res = tgates.CCCRy(res.states[-1], 0, 1, 3, 2, thetay00)
    res = tgates.X(res.states[-1], 3)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, thetay11)

    res = tgates.X(res.states[-1], 1)
    return tgates.X(res.states[-1], 0)

def Kb1d(psi0):
    thetay00 = 1.85806
    thetay10 = 2.48274
    thetay11 = 1.5708

    res = tgates.X(psi0, 0)
    res = tgates.X(res.states[-1], 1)

    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, -thetay11)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, -thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.X(res.states[-1], 3)
    res = tgates.CCCRy(res.states[-1], 0, 1, 3, 2, -thetay00)
    res = tgates.X(res.states[-1], 3)

    res = tgates.X(res.states[-1], 1)
    return tgates.X(res.states[-1], 0)

def Kb2(psi0):
    # \[Theta]y00 -> 0.554811, \[Theta]y10 -> \ 0.405465, \[Theta]y11 -> 1.5708
    thetay00 = 0.554811
    thetay10 = 0.405465
    thetay11 = 1.5708

    res = tgates.X(psi0, 3)
    res = tgates.CRy(res.states[-1], 3, 2, thetay00)
    res = tgates.X(res.states[-1], 3)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CRy(res.states[-1], 2, 3, thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.CRy(res.states[-1], 2, 3, thetay11)


    res = tgates.X(res.states[-1], 0)
    res = tgates.X(res.states[-1], 1)

    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, -thetay11)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, -thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.X(res.states[-1], 3)
    res = tgates.CCCRy(res.states[-1], 0, 1, 3, 2, -thetay00)
    res = tgates.X(res.states[-1], 3)

    res = tgates.X(res.states[-1], 1)
    return tgates.X(res.states[-1], 0)

def Kb2d(psi0):
    thetay00 = 0.554811
    thetay10 = 0.405465
    thetay11 = 1.5708

    res = tgates.X(psi0, 0)
    res = tgates.X(res.states[-1], 1)

    res = tgates.X(res.states[-1], 3)
    res = tgates.CCCRy(res.states[-1], 0, 1, 3, 2, thetay00)
    res = tgates.X(res.states[-1], 3)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.CCCRy(res.states[-1], 0, 1, 2, 3, thetay11)

    res = tgates.X(res.states[-1], 1)
    res = tgates.X(res.states[-1], 0)


    res = tgates.CRy(res.states[-1], 2, 3, -thetay11)

    res = tgates.X(res.states[-1], 2)
    res = tgates.CRy(res.states[-1], 2, 3, -thetay10)
    res = tgates.X(res.states[-1], 2)

    res = tgates.X(res.states[-1], 3)
    res = tgates.CRy(res.states[-1], 3, 2, -thetay00)
    return tgates.X(res.states[-1], 3)

def Dd(psi0):
    return tgates.CP(psi0, 2, 3, np.pi, b=0b00)

def reg_SWAP(psi0):
    res = tgates.SWAP(psi0, 1, 3)
    return tgates.SWAP(res.states[-1], 0, 2)


'''
# El algoritmo
# Estado fiducial
print('Preparando estado fiducial...\t {}/{}/{} - {}:{}:{}'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
psi0 = tensor(basis(2,0), basis(2,0), basis(2,0), basis(2,0))

# Preparación del estado inicial
print('Preparando estado inicial (1/4)...\t {}/{}/{} - {}:{}:{}'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = tgates.H(psi0,0)
res = tgates.H(res.states[-1],1)
print('Preparando estado inicial (2/4)...\t {}/{}/{} - {}:{}:{}'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = Kb1(res.states[-1])
print('Preparando estado inicial (3/4)...\t {}/{}/{} - {}:{}:{}'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = Kb2(res.states[-1])
'''


res = qload('itj_33')

for i in range(33,60):
    
    print('Iteración {}: Aplicando Kb2^dag...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb2d(res.states[-1])
    print('Iteración {}: Aplicando Kb1^dag...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb1d(res.states[-1])
    print('Iteración {}: Aplicando Dd...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Dd(res.states[-1])
    print('Iteración {}: Aplicando Kb1...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb1(res.states[-1])
    print('Iteración {}: Aplicando Kb2...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb2(res.states[-1])
    print('Iteración {}: Aplicando reg_SWAP...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = reg_SWAP(res.states[-1])
    
    print('Iteración {}: Guardando resultado de iteración...\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    qsave(res, 'itj_{}'.format(i+1))

    print('Iteración {}: Terminada\t {}/{}/{} - {}:{}:{}'.format(i+1, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))


''' # +8
for i in range(21):
    
    print('Iteración {}: Aplicando T1...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = T1(res.states[-1])
    print('Iteración {}: Aplicando Kb2^dag...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb2d(res.states[-1])
    print('Iteración {}: Aplicando Kb1^dag...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb1d(res.states[-1])
    print('Iteración {}: Aplicando Dd...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Dd(res.states[-1])
    print('Iteración {}: Aplicando Kb1...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb1(res.states[-1])
    print('Iteración {}: Aplicando Kb2...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = Kb2(res.states[-1])
    print('Iteración {}: Aplicando T1^dag...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = T1d(res.states[-1])
    print('Iteración {}: Aplicando reg_SWAP...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    res = reg_SWAP(res.states[-1])
    
    print('Iteración {}: Guardando resultado de iteración...\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
    qsave(res, 'itj_{}'.format(i+8))

    print('Iteración {}: Terminada\t {}/{}/{} - {}:{}:{}'.format(i+8, time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))

'''
