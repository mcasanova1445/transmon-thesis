# Importando librerias y modulos necesarios

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *

class QuantumChip:
    def __init__(self,
                 N = 50,
                 wq = [5.0 * 2 * np.pi, 6.0 * 2 * np.pi, 7.0 * 2 * np.pi, 8.0 * 2 * np.pi],
                 wr = 10.0 * 2 * np.pi,
                 wq_swap = 9 * 2 * np.pi,
                 g = [0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi],
                 kappa = 0.001,
                 gamma = [5e-6, 5e-6, 5e-6, 5e-6],
                 psi0 = [[1,0], [1,0], [1,0], [1,0]]):

        if(len(wq) == len(g) and
           (len(wq) == len(gamma) or len(gamma) == 0) and
           len(wq) == len(psi0) and
           all(len(psi0[partition]) == 2 for partition in range(len(psi0)))):

            # Parametros del sistema
            self.N = N
            self.qN = len(wq)

            self.wr = wr
            self.wq = np.array(wq)
            self.wq_swap = wq_swap

            self.g = np.array(g)

            self.D = self.wq - self.wr
            self.D_swap = self.wq_swap - self.wr

            self.chi = self.g**2 / abs(self.wr-self.wq)

            self.kappa = kappa
            self.gamma = np.array(gamma)
            self.c_ops = []

            for i in range(len(gamma)):
                self.c_ops = self.c_ops + [np.sqrt(self.gamma[i]) * self.qop('sm', i)]

            self.psi0 = (psi0[0][0]*basis(2, 0) + psi0[0][1]*basis(2, 1)).unit()

            for i in range(self.qN-1):
                self.psi0 = tensor(self.psi0, (psi0[i+1][0]*basis(2, 0) + psi0[i+1][1]*basis(2, 1)).unit() )
                
            '''
            # Operadores del resonador
            self.a = destroy(self.N)
            self.n = a.dag() * a
            self.Id_r = qeye(self.N)
            '''

        else:
            print('Dimensiones inconsistentes')

    # Funcion para graficar
    # TODO elegir colores lindos para las gráficas

    def plot_drive_expect(self,res,args):
        tlist = res.times

        if args == 0:
            fig, axes = plt.subplots(1, 1, sharex=True, figsize=(12,4))

            for i in range(self.qN):
                axes.plot(tlist, np.real(expect(self.qop('n',i), res.states)), '{}'.format(i/self.qN), linewidth=2, label="qubit {}".format(i))
            axes.set_ylim(0, 1)

            axes.set_xlabel("Time (ns)", fontsize=16)
            axes.set_ylabel("Occupation probability", fontsize=16)
            axes.legend()

        else:
            fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12,8))

            axes[0].plot(tlist, np.array(list(self.ksi_t(tlist,args))) / (2*np.pi), 'b', linewidth=2, label="drive envelope")
            axes[0].set_ylabel("Energy (GHz)", fontsize=16)
            axes[0].legend()

            for i in range(self.qN):
                axes.plot(tlist, np.real(expect(self.qop('n',i), res.states)), '{}'.format(i/self.qN), linewidth=2, label="qubit {}".format(i))
            axes[1].set_ylim(0, 1)

            axes[1].set_xlabel("Time (ns)", fontsize=16)
            axes[1].set_ylabel("Occupation probability", fontsize=16)
            axes[1].legend()

        fig.tight_layout()


    # Funciones para generar los operadores de los transmones

    def qop_part(self, operator, target):
        if target == 0:
            qop_dict = {'sm' : destroy(2), 'sp' : (destroy(2)).dag(), 
                        'sx' : sigmax(), 'sy' : sigmay(), 'sz' : sigmaz(),
                        'n' : (destroy(2)).dag() * destroy(2)}
            qop_part_out = qop_dict[operator]
        else:
            qop_part_out = qeye(2)

        return qop_part_out

    def qop(self, operator, target):
        qop_out = self.qop_part(operator, target)

        for i in range(self.qN-1):
            qop_out = tensor(qop_out, self.qop_part(operator, target-1-i))

        return qop_out
        
    # Compuertas

    def Rx(self, target, theta):
        tlist = np.linspace(0, 10, 200)

        wd = self.wq[target]

        Dr = self.wr-wd
        Dq = self.wq-wd

        Hsyst = 0
        for i in range(self.qN):
            Hsyst = Hsyst - Dq[i]*self.qop('sz',i)/2

        H_t = [[self.qop('sx',target)/2, ksi_t], Hsyst]

        args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : self.wq[target]}
        res = mesolve(H_t, self.psi0, tlist, self.c_ops, [], args = args)
        self.psi0 = res.states[-1]

        # self.plot_drive_expect(res,args)

        return res

    def Ry(self, target, theta):
        tlist = np.linspace(0, 10, 200)

        wd = self.wq[target]

        Dr = self.wr-wd
        Dq = self.wq-wd

        Hsyst = 0
        for i in range(self.qN):
            Hsyst = Hsyst - Dq[i]*self.qop('sz',i)/2

        H_t = [[self.qop('sy',target)/2, ksi_t], Hsyst]

        args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : self.wq[target]}
        res = mesolve(H_t, self.psi0, tlist, self.c_ops, [], args = args)
        self.psi0 = res.states[-1]

        # self.plot_drive_expect(res,args)

        return res

    def Rz(self, target, theta):
        res = self.Ry(target, np.pi/2)
        res = self.Rx(target, theta)
        return self.Ry(target, -np.pi/2)

    def X(self, target):
        return self.Rx(target, np.pi)

    def Y(self, target):
        return self.Ry(target, np.pi)

    def Z(self, target, theta):
        return self.Rz(target, np.pi)

    def H(self, target):
        res = self.Ry(target, np.pi/2)
        return self.X(target)

    def sqrtiSWAP(self, target1, target2):
        wq = self.wq
        
        wq[target1] = self.wq_swap
        wq[target2] = self.wq_swap

        D = wq - self.wr

        J = np.abs(self.g[target1] * self.g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

        tf = np.pi/(4*J)
        tlist = np.linspace(0, tf, 250)

        Hsyst = self.g[target1]*self.g[target2] * (self.qop('sp',target1)*self.qop('sm',target2) + self.qop('sm',target1)*self.qop('sp',target2)) / (self.D_swap)

        res = mesolve(Hsyst, self.psi0, tlist, self.c_ops, [])
        self.psi0 = res.states[-1]

        args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]}
        # self.plot_drive_expect(res,args)

        return res

    def iSWAP(self, target1, target2):
        wq = self.wq
        
        wq[target1] = self.wq_swap
        wq[target2] = self.wq_swap

        D = wq - self.wr

        J = np.abs(self.g[target1] * self.g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

        tf = np.pi/(2*J)
        tlist = np.linspace(0, tf, 500)

        Hsyst = self.g[target1]*self.g[target2] * (self.qop('sp',target1)*self.qop('sm',target2) + self.qop('sm',target1)*self.qop('sp',target2)) / (self.D_swap)

        res = mesolve(Hsyst, self.psi0, tlist, self.c_ops, [])
        self.psi0 = res.states[-1]

        args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]} # Revisar este wq si decido graficar otra vez
        # self.plot_drive_expect(res,args)

        return res

    def CNOT(self, control, target):
        res = self.H(target)
        res = self.Rz(target, -np.pi/2)
        res = self.Rz(control, -np.pi/2)
        res = self.iSWAP(control, target)
        res = self.H(control)
        res = self.iSWAP(control, target)
        res = self.Rx(target, np.pi/2)
        res = self.iSWAP(control, target)
        res = self.Rx(control, np.pi/2)
        res = self.iSWAP(control, target)
        return self.Rx(target, np.pi/2)

    def CRy(self, control, target, theta):
        res = self.Ry(target, theta/2)
        res = self.CNOT(control, target)
        res = self.Ry(target, -theta/2)
        return self.CNOT(control, target)

    def CRz(self, control, target, theta):
        res = self.Rz(target, theta/2)
        res = self.CNOT(control, target)
        res = self.Rz(target, -theta/2)
        return self.CNOT(control, target)

    def SWAP(self, target1, target2):
        res = self.CNOT(target1, target2)
        res = self.CNOT(target2, target1)
        return self.CNOT(target1, target2)

    def CH(self, control, target):
        res = self.Ry(target, np.pi/4)
        res = self.CNOT(control, target)
        return self.Ry(target, -np.pi/4)

    def CP(self, control, target, theta, b = 0b11):
        if b == 0b00:
            res = self.Rz(control, -3*theta/4)
            res = self.Rz(target, -3*theta/4)
            res = self.CRz(control, target, theta/2)
            res = self.CRz(target, control, theta/2)

        elif b == 0b01:
            res = self.Rz(control, -3*theta/4)
            res = self.Rz(target, 5*theta/4)
            res = self.CRz(control, target, -3*theta/2)
            res = self.CRz(target, control, theta/2)

        elif b == 0b10:
            res = self.Rz(control, theta/4)
            res = self.Rz(target, theta/4)
            res = self.CRz(control, target, -3*theta/2)
            res = self.CRz(target, control, theta/2)

        elif b == 0b11:
            res = self.Rz(control, theta/4)
            res = self.Rz(target, theta/4)
            res = self.CRz(control, target, theta/2)
            res = self.CRz(target, control, theta/2)

        return res

    def Toffoli(self, control1, control2, target):
        res = self.H(target)
        res = self.CRz(control2, target, -np.pi/2)
        res = self.CNOT(control1, control2)
        res = self.CRz(control2, target, np.pi/2)
        res = self.CNOT(control1, control2)
        res = self.CRz(control1, target, -np.pi/2)
        res = self.H(target)
        return self.CP(control1, control2, -np.pi/2, b = 0b11)

    def CCRz(self, control1, control2, target, theta):
        res = self.CRz(control2, target, theta/2)
        res = self.CNOT(control1, control2)
        res = self.CRz(control2, target, -theta/2)
        res = self.CNOT(control1, control2)
        return self.CRz(control1, target, theta/2)

    def CCRy(self, control1, control2, target, theta):
        res = self.CRy(control2, target, theta/2)
        res = self.CNOT(control1, control2)
        res = self.CRy(control2, target, -theta/2)
        res = self.CNOT(control1, control2)
        return self.CRy(control1, target, theta/2)

    def CCP(self, control1, control2, target, theta, b = 0b11):
        res = self.CP(control2, target, theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
        res = self.CNOT(control1, control2)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
        res = self.CP(control2, target, -theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
        res = self.CNOT(control1, control2)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
        return self.CP(control1, target, theta/2, b = b)

    def CCNOT(self, control1, control2, target):
        return self.Toffoli(control1, control2, target)

    def mZ(self, target):
        res = self.Ry(target, np.pi)
        return self.Rx(target, np.pi)

    def CCCNOT(self, control1, control2, control3, target):
        res = self.H(target)
        res = self.CCRz(control2, control3, target, -np.pi/2)
        res = self.CNOT(control1, control2)
        res = self.CCRz(control2, control3, target, np.pi/2)
        res = self.CNOT(control1, control2)
        res = self.CCRz(control1, control3, target, -np.pi/2)
        res = self.H(target)
        res = self.CP(control2, control3, -np.pi/4)
        res = self.CNOT(control1, control2)
        res = self.CP(control2, control3, np.pi/4)
        res = self.CNOT(control1, control2)
        return self.CP(control1, control3, -np.pi/4)

    def CCCRy(self, control1, control2, control3, target, theta):
        res = self.CRy(control3, target, theta/2)
        res = self.CCNOT(control1, control2, control3)
        res = self.CRy(control3, target, -theta/2)
        res = self.CCNOT(control1, control2, control3)
        return self.CCRy(control1, control2, target, theta/2)

    def CCCRz(self, control1, control2, control3, target, theta):
        res = self.CRz(control3, target, theta/2)
        res = self.CCNOT(control1, control2, control3)
        res = self.CRz(control3, target, -theta/2)
        res = self.CCNOT(control1, control2, control3)
        return self.CCRz(control1, control2, target, theta/2)

    def CCCP(self, control1, control2, control3, target, theta, b = 0b11):
        res = self.CP(control3, target, theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
            res = self.X(control2)
        res = self.CCNOT(control1, control2, control3)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
            res = self.X(control2)
        res = self.CP(control3, target, -theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
            res = self.X(control2)
        res = self.CCNOT(control1, control2, control3)
        if b == 0b00 or b == 0b01:
            res = self.X(control1)
            res = self.X(control2)
        return self.CCP(control1, control2, target, theta/2, b = b)


# Funciones de los pulsos

def gaussianpulse(x,ts,tf):
    s = (tf-ts)/6
    m = (ts+tf)/2
    return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))*norm.pdf(x, loc = m, scale = s)/0.997300204

def squarepulse(x,ts,tf):
    s = (tf-ts)/6
    m = (ts+tf)/2
    return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))/(6*s)

# Funciones para la dependencia temporal de los Hamiltonianos

def ksi_t(t, args):
    return args['A'] * gaussianpulse(t,args['ts'],args['tf'])

def ksi_tm(t, args):
    return args['A'] * gaussianpulse(t,args['ts'],args['tf']) * np.exp(-1j*args['w']*(t-args['ts']))

def ksi_tp(t, args):
    return args['A'] * gaussianpulse(t,args['ts'],args['tf']) * np.exp(1j*args['w']*(t-args['ts']))

def ksiS_t(t, args):
    return args['A'] * squarepulse(t,args['ts'],args['tf'])

def ksiS_tm(t, args):
    return args['A'] * squarepulse(t,args['ts'],args['tf']) * np.exp(-1j*args['w']*(t-args['ts']))

def ksiS_tp(t, args):
    return args['A'] * squarepulse(t,args['ts'],args['tf']) * np.exp(1j*args['w']*(t-args['ts']))


