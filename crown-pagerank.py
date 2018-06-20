# Importando librerías necesarias
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from qutip import *
import tgates

# Estado inicial
psi0 = tensor(basis(2,0),
              basis(2,0),
              basis(2,0),
              basis(2,0))

# Definiendo operadores
def R2(psi, b0, b1):
    res = tgates.CNOT(psi, b1, b0)
    return tgates.X(res.states[-1], b1)



# El algoritmo


