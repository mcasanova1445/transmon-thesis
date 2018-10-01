import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *
import tgates8
import time

def QFT8(psi0):
    res = tgates8.H(psi0, 0)
    res = tgates8.CP(res.states[-1], 1, 0, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 2, 0, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], 3, 0, 2*np.pi/2**4)
    res = tgates8.CP(res.states[-1], 4, 0, 2*np.pi/2**5)
    res = tgates8.CP(res.states[-1], 5, 0, 2*np.pi/2**6)
    res = tgates8.CP(res.states[-1], 6, 0, 2*np.pi/2**7)
    res = tgates8.CP(res.states[-1], 7, 0, 2*np.pi/2**8)
    res = tgates8.H(res.states[-1], 1)
    res = tgates8.CP(res.states[-1], 2, 1, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 3, 1, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], 4, 1, 2*np.pi/2**4)
    res = tgates8.CP(res.states[-1], 5, 1, 2*np.pi/2**5)
    res = tgates8.CP(res.states[-1], 6, 1, 2*np.pi/2**6)
    res = tgates8.CP(res.states[-1], 7, 1, 2*np.pi/2**7)
    res = tgates8.H(res.states[-1], 2)
    res = tgates8.CP(res.states[-1], 3, 2, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 4, 2, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], 5, 2, 2*np.pi/2**4)
    res = tgates8.CP(res.states[-1], 6, 2, 2*np.pi/2**5)
    res = tgates8.CP(res.states[-1], 7, 2, 2*np.pi/2**6)
    res = tgates8.H(res.states[-1], 3)
    res = tgates8.CP(res.states[-1], 4, 3, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 5, 3, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], 6, 3, 2*np.pi/2**4)
    res = tgates8.CP(res.states[-1], 7, 3, 2*np.pi/2**5)
    res = tgates8.H(res.states[-1], 4)
    res = tgates8.CP(res.states[-1], 5, 4, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 6, 4, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], 7, 4, 2*np.pi/2**4)
    res = tgates8.H(res.states[-1], 5)
    res = tgates8.CP(res.states[-1], 6, 5, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], 7, 5, 2*np.pi/2**3)
    res = tgates8.H(res.states[-1], 6)
    res = tgates8.CP(res.states[-1], 7, 6, 2*np.pi/2**2)
    return tgates8.H(res.states[-1], 7)

def QFT4(psi0, target1, target2, target3, target4):
    res = tgates8.H(psi0, target1)
    res = tgates8.CP(res.states[-1], target2, target1, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], target3, target1, 2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], target4, target1, 2*np.pi/2**4)
    res = tgates8.H(res.states[-1], target2)
    res = tgates8.CP(res.states[-1], target3, target2, 2*np.pi/2**2)
    res = tgates8.CP(res.states[-1], target4, target2, 2*np.pi/2**3)
    res = tgates8.H(res.states[-1], target3)
    res = tgates8.CP(res.states[-1], target4, target3, 2*np.pi/2**2)
    res = tgates8.H(res.states[-1], target4)
    res = tgates8.SWAP(res.states[-1], target1, target4)
    return tgates8.SWAP(res.states[-1], target2, target3)

def QFT4d(psi0, target1, target2, target3, target4):
    res = tgates8.SWAP(psi0, target2, target3)
    res = tgates8.SWAP(res.states[-1], target1, target4)
    res = tgates8.H(res.states[-1], target4)
    res = tgates8.CP(res.states[-1], target4, target3, -2*np.pi/2**2)
    res = tgates8.H(res.states[-1], target3)
    res = tgates8.CP(res.states[-1], target4, target2, -2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], target3, target2, -2*np.pi/2**2)
    res = tgates8.H(res.states[-1], target2)
    res = tgates8.CP(res.states[-1], target4, target1, -2*np.pi/2**4)
    res = tgates8.CP(res.states[-1], target3, target1, -2*np.pi/2**3)
    res = tgates8.CP(res.states[-1], target2, target1, -2*np.pi/2**2)
    return tgates8.H(res.states[-1], target1)

def CSWAP(psi0, control, target1, target2):
    res = tgates8.CCNOT(psi0, control, target1, target2)
    res = tgates8.CCNOT(psi0, control, target2, target1)
    return tgates8.CCNOT(psi0, control, target1, target2)

def MUL7(psi0, target1, target2, target3, target4):
    res = tgates8.X(psi0, target1)
    res = tgates8.X(res.states[-1], target2)
    res = tgates8.X(res.states[-1], target3)
    res = tgates8.X(res.states[-1], target4)
    res = tgates8.SWAP(res.states[-1], target2, target3)
    res = tgates8.SWAP(res.states[-1], target1, target2)
    return tgates8.SWAP(res.states[-1], target1, target4)

def CMUL7(psi0, control, target1, target2, target3, target4):
    res = tgates8.CNOT(psi0, control, target1)
    res = tgates8.CNOT(res.states[-1], control, target2)
    res = tgates8.CNOT(res.states[-1], control, target3)
    res = tgates8.CNOT(res.states[-1], control, target4)
    res = CSWAP(res.states[-1], control, target2, target3)
    res = CSWAP(res.states[-1], control, target1, target2)
    return CSWAP(res.states[-1], control, target1, target4)

def MUL3(psi0, target1, target2, target3, target4):
    res = tgates8.X(psi0, target1)
    res = tgates8.CNOT(res.states[-1], target1, target4)
    res = tgates8.CCCNOT(res.states[-1], target1, target4, target3, target2)
    res = tgates8.CNOT(res.states[-1], target1, target4)
    res = tgates8.CCNOT(res.states[-1], target1, target4, target3)
    return tgates8.X(res.states[-1], target1)

def CMUL3(psi0, target1, target2, target3, target4):
    res = tgates8.CNOT(psi0, control, control, target1)
    res = tgates8.CCNOT(res.states[-1], control, target1, target4)
    res = tgates8.CCCCNOT(res.states[-1], control, target1, target4, target3, target2)
    res = tgates8.CCNOT(res.states[-1], control, target1, target4)
    res = tgates8.CCCNOT(res.states[-1], control, target1, target4, target3)
    return tgates8.CNOT(res.states[-1], control, target1)

qN = 2**8


# El algoritmo
# Estado fiducial
print('{}/{}/{} - {}:{}:{}\t Preparando estado fiducial...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
psi0 = tensor(basis(2,0), basis(2,0), basis(2,0), basis(2,0), basis(2,0), basis(2,0), basis(2,0), basis(2,0))

# Preparación del estado inicial
print('{}/{}/{} - {}:{}:{}\t Preparando estado inicial...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = tgates8.X(psi0, 7)

qsave(res, 'rp_0')

# Transformada de Hadamard en el primer registro
print('{}/{}/{} - {}:{}:{}\t Aplicando transformada de Hadamard al primer registro...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = tgates8.H(res.states[-1], 0)
res = tgates8.H(res.states[-1], 1)
res = tgates8.H(res.states[-1], 2)
res = tgates8.H(res.states[-1], 3)

qsave(res, 'rp_1')

# Aplicando CMUL3(c=3)^1
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=3)^1 (1/1)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 3, 4, 5, 6, 7)

qsave(res, 'rp_2')

# Aplicando CMUL3(c=2)^2
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=2)^2 (1/2)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 2, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=2)^2 (2/2)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 2, 4, 5, 6, 7)

qsave(res, 'rp_3')

# Aplicando CMUL3(c=1)^4
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=1)^4 (1/4)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 1, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=1)^4 (2/4)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 1, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=1)^4 (3/4)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 1, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=1)^4 (4/4)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 1, 4, 5, 6, 7)

qsave(res, 'rp_4')

# Aplicando CMUL3(c=0)^8
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (1/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (2/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (3/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (4/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (5/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (6/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (7/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)
print('{}/{}/{} - {}:{}:{}\t Aplicando CMUL3(c=0)^8 (8/8)...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = CMUL3(res.states[-1], 0, 4, 5, 6, 7)

qsave(res, 'rp_5')

# Aplicando la transformada cuántica inversa de Fourier sobre el primer registro
print('{}/{}/{} - {}:{}:{}\t Aplicando QFT4d...'.format(time.localtime()[0], time.localtime()[1], time.localtime()[2], time.localtime()[3], time.localtime()[4], time.localtime()[5]))
res = QFT4d(res.states[-1], 0, 1, 2, 3)

qsave(res, 'rp_6')
