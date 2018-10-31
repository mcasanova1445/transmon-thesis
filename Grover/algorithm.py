# PageRank del grafo corona

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *
import tgates
import time

def Htrans(psi0):
    res = tgates.H(psi0, 0)
    res = tgates.H(res.states[-1], 1)
    res = tgates.H(res.states[-1], 2)
    return tgates.H(res.states[-1], 3)

def Us(psi0):
    res = Htrans(psi0)
    res = tgates.CCCP(res.states[-1], 1, 2, 3, 0, np.pi, b = 0b00)
    return Htrans(res.states[-1])

def Uomega(psi0):
    res = tgates.X(psi0, 0)
    res = tgates.CCCP(res.states[-1], 0, 1, 2, 3, np.pi, b = 0b10)
    return tgates.X(res.states[-1], 0)

'''
def Uomega(psi0):
    return tgates.CCP(psi0, 0, 3, 2, np.pi, b = 0b10)
'''

'''
def Uomega(psi0):
    return tgates.CP(psi0, 2, 1, np.pi, b = 0b01)
'''

qN = 2**4


# El algoritmo
# Estado fiducial
print('{}/{}/{} - {}:{}:{}\t Preparando estado fiducial...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
psi0 = tensor(basis(2,0), basis(2,0), basis(2,0), basis(2,0))

# Preparación del estado inicial
print('{}/{}/{} - {}:{}:{}\t Preparando estado inicial...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = Htrans(psi0)

Nit = 2*int(np.pi*np.sqrt(qN)/4)+1

for i in range(Nit):
    
    print('{}/{}/{} - {}:{}:{}\t Iteración {}/{}: Aplicando Uomega...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5], i+1, Nit))
    res = Uomega(res.states[-1])
    print('{}/{}/{} - {}:{}:{}\t Iteración {}/{}: Aplicando Us...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5], i+1, Nit))
    res = Us(res.states[-1])
    
    print('{}/{}/{} - {}:{}:{}\t Iteración {}/{}: Guardando resultado de iteración...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5], i+1, Nit))
    qsave(res, 'it_{}'.format(i+1))

    print('{}/{}/{} - {}:{}:{}\t Iteración {}/{}: Terminada'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5], i+1, Nit))

