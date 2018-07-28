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
                 psi0 = tensor(basis(2,0), basis(2,0), basis(2,0), basis(2,0)) ):

        if(len(wq) == len(g) and
           (len(wq) == len(gamma) or len(gamma) == 0) and
           len(wq) == len(psi0)):

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
                self.c_ops = self.c_ops + [np.sqrt(self.gamma[i]) * qop('sm', i)]

            self.psi0 = basis(2, psi0[0])

            for i in range(self.qN-1):
                self.psi0 = tensor(self.psi, basis(2, psi0[i+1]))
        
        else:
            print('Dimensiones inconsistentes')

    # Funciones de los pulsos

    def gaussianpulse(x,ts,tf):
        s = (tf-ts)/6
        m = (ts+tf)/2
        return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))*norm.pdf(x, loc = m, scale = s)/0.997300204

    def squarepulse(x,ts,tf):
        s = (tf-ts)/6
        m = (ts+tf)/2
        return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))/(6*s)

    # Funcion para graficar
    # TODO elegir colores lindos para las gráficas

    def plot_drive_expect(res,args):
        tlist = res.times

        if args == 0:
            fig, axes = plt.subplots(1, 1, sharex=True, figsize=(12,4))

            for i in range(self.qN):
                axes.plot(tlist, np.real(expect(qop('n',i), res.states)), '{}'.format(i/self.qN), linewidth=2, label="qubit {}".format(i))
            axes.set_ylim(0, 1)

            axes.set_xlabel("Time (ns)", fontsize=16)
            axes.set_ylabel("Occupation probability", fontsize=16)
            axes.legend()

        else:
            fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12,8))

            axes[0].plot(tlist, np.array(list(ksi_t(tlist,args))) / (2*np.pi), 'b', linewidth=2, label="drive envelope")
            axes[0].set_ylabel("Energy (GHz)", fontsize=16)
            axes[0].legend()

            for i in range(self.qN):
                axes.plot(tlist, np.real(expect(qop('n',i), res.states)), '{}'.format(i/self.qN), linewidth=2, label="qubit {}".format(i))
            axes[1].set_ylim(0, 1)

            axes[1].set_xlabel("Time (ns)", fontsize=16)
            axes[1].set_ylabel("Occupation probability", fontsize=16)
            axes[1].legend()

        fig.tight_layout()


    # Operadores del resonador

    a = destroy(self.N)
    n = a.dag() * a
    Id_r = qeye(self.N)

    # Funciones para generar los operadores de los transmones

    def qop_part(operator, target):
        if target == 0:
            qop_dict = {'sm' : destroy(2), 'sp' : (destroy(2)).dag(), 
                        'sx' : sigmax(), 'sy' : sigmay(), 'sz' : sigmaz(),
                        'n' : (destroy(2)).dag() * destroy(2)}
            qop_part_out = qop_dict[operator]
        else:
            qop_part_out = qeye(2)

        return qop_part_out

    def qop(operator, target):
        qop_out = qop_part(operator, target)

        for i in range(self.qN-1):
            qop_out = tensor(qop_out, qop_part(operator, target-1-i))

        return qop_out
        
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
        return args['A'] * np.exp(-1j*args['w']*(t-args['ts']))

    def ksiS_tp(t, args):
        return args['A'] * np.exp(1j*args['w']*(t-args['ts']))

    # Compuertas

    def Rx(self, target, theta):
        tlist = np.linspace(0, 10, 200)

        wd = self.wq[target]

        Dr = self.wr-wd
        Dq = self.wq-wd

        Hsyst = 0
        for i in range(self.qN):
            Hsyst = Hsyst - Dq[i]*qop('sz',i)/2

        H_t = [[qop('sx',target)/2, ksi_t], Hsyst]

        args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : self.wq[target]}
        res = mesolve(H_t, self.psi0, tlist, self.c_ops, [], args = args)
        self.psi0 = res.states[-1]

        # plot_drive_expect(res,args)

        return res

    def Ry(self, target, theta):
        tlist = np.linspace(0, 10, 200)

        wd = self.wq[target]

        Dr = self.wr-wd
        Dq = self.wq-wd

        Hsyst = 0
        for i in range(self.qN):
            Hsyst = Hsyst - Dq[i]*qop('sz',i)/2

        H_t = [[qop('sy',target)/2, ksi_t], Hsyst]

        args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : self.wq[target]}
        res = mesolve(H_t, self.psi0, tlist, self.c_ops, [], args = args)
        self.psi0 = res.states[-1]

        # plot_drive_expect(res,args)

        return res

    def Rz(self, target, theta):
        res = Ry(self.psi0, target, np.pi/2)
        res = Rx(self.psi0, target, theta)
        return Ry(self.psi0, target, -np.pi/2)

    def X(self, target):
        return Rx(self.psi0, target, np.pi)

    def Y(self, target):
        return Ry(self.psi0, target, np.pi)

    def Z(self, target, theta):
        return Rz(self.psi0, target, np.pi)

    def H(self, target):
        res = Ry(self.psi0, target, np.pi/2)
        return X(self.psi0, target)

    def sqrtiSWAP(self, target1, target2):
        wq = self.wq
        
        wq[target1] = self.wq_swap
        wq[target2] = self.wq_swap

        D = wq - self.wr

        J = np.abs(self.g[target1] * self.g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

        tf = np.pi/(4*J)
        tlist = np.linspace(0, tf, 250)

        Hsyst = self.g[target1]*self.g[target2] * (qop('sp',target1)*qop('sm',target2) + qop('sm',target1)*qop('sp',target2)) / (self.D_swap)

        res = mesolve(Hsyst, self.psi0, tlist, self.c_ops, [])
        self.psi0 = res.states[-1]

        args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]}
        # plot_drive_expect(res,args)

        return res

    def iSWAP(self, target1, target2):
        wq = self.wq
        
        wq[target1] = self.wq_swap
        wq[target2] = self.wq_swap

        D = wq - self.wr

        J = np.abs(self.g[target1] * self.g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

        tf = np.pi/(2*J)
        tlist = np.linspace(0, tf, 500)

        Hsyst = self.g[target1]*self.g[target2] * (qop('sp',target1)*qop('sm',target2) + qop('sm',target1)*qop('sp',target2)) / (self.D_swap)

        res = mesolve(Hsyst, self.psi0, tlist, self.c_ops, [])
        self.psi0 = res.states[-1]

        args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]} # Revisar este wq si decido graficar otra vez
        # plot_drive_expect(res,args)

        return res

    def CNOT(self, control, target):
        res = H(self.psi0, target)
        res = Rz(self.psi0, target, -np.pi/2)
        res = Rz(self.psi0, control, -np.pi/2)
        res = iSWAP(self.psi0, control, target)
        res = H(self.psi0, control)
        res = iSWAP(self.psi0, control, target)
        res = Rx(self.psi0, target, np.pi/2)
        res = iSWAP(self.psi0, control, target)
        res = Rx(self.psi0, control, np.pi/2)
        res = iSWAP(self.psi0, control, target)
        return Rx(self.psi0, target, np.pi/2)

    def CRy(self, control, target, theta):
        res = Ry(self.psi0,target,theta/2)
        res = CNOT(self.psi0,control,target)
        res = Ry(self.psi0,target,-theta/2)
        return CNOT(self.psi0,control,target)

    def CRz(self, control, target, theta):
        res = Rz(self.psi0,target,theta/2)
        res = CNOT(self.psi0,control,target)
        res = Rz(self.psi0,target,-theta/2)
        return CNOT(self.psi0,control,target)

    def SWAP(self, target1, target2):
        res = CNOT(self.psi0, target1, target2)
        res = CNOT(self.psi0, target2, target1)
        return CNOT(self.psi0, target1, target2)

    def CH(self, control, target):
        res = Ry(self.psi0, target, np.pi/4)
        res = CNOT(self.psi0, control, target)
        return Ry(self.psi0, target, -np.pi/4)

    def CP(self, control, target, theta, b = 0b11):
        if b == 0b00:
            res = Rz(self.psi0, control, -3*theta/4)
            res = Rz(self.psi0, target, -3*theta/4)
            res = CRz(self.psi0, control, target, theta/2)
            res = CRz(self.psi0, target, control, theta/2)

        elif b == 0b01:
            res = Rz(self.psi0, control, -3*theta/4)
            res = Rz(self.psi0, target, 5*theta/4)
            res = CRz(self.psi0, control, target, -3*theta/2)
            res = CRz(self.psi0, target, control, theta/2)

        elif b == 0b10:
            res = Rz(self.psi0, control, theta/4)
            res = Rz(self.psi0, target, theta/4)
            res = CRz(self.psi0, control, target, -3*theta/2)
            res = CRz(self.psi0, target, control, theta/2)

        elif b == 0b11:
            res = Rz(self.psi0, control, theta/4)
            res = Rz(self.psi0, target, theta/4)
            res = CRz(self.psi0, control, target, theta/2)
            res = CRz(self.psi0, target, control, theta/2)

        return res

    def Toffoli(self, control1, control2, target):
        res = H(self.psi0, target)
        res = CRz(self.psi0, control2, target, -np.pi/2)
        res = CNOT(self.psi0, control1, control2)
        res = CRz(self.psi0, control2, target, np.pi/2)
        res = CNOT(self.psi0, control1, control2)
        res = CRz(self.psi0, control1, target, -np.pi/2)
        res = H(self.psi0, target)
        return CP(self.psi0, control1, control2, -np.pi/2, b = 0b11)

    def CCRz(self, control1, control2, target, theta):
        res = CRz(self.psi0, control2, target, theta/2)
        res = CNOT(self.psi0, control1, control2)
        res = CRz(self.psi0, control2, target, -theta/2)
        res = CNOT(self.psi0, control1, control2)
        return CRz(self.psi0, control1, target, theta/2)

    def CCRy(self, control1, control2, target, theta):
        res = CRy(self.psi0, control2, target, theta/2)
        res = CNOT(self.psi0, control1, control2)
        res = CRy(self.psi0, control2, target, -theta/2)
        res = CNOT(self.psi0, control1, control2)
        return CRy(self.psi0, control1, target, theta/2)

    def CCP(self, control1, control2, target, theta, b = 0b11):
        res = CP(self.psi0, control2, target, theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
        res = CNOT(self.psi0, control1, control2)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
        res = CP(self.psi0, control2, target, -theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
        res = CNOT(self.psi0, control1, control2)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
        return CP(self.psi0, control1, target, theta/2, b = b)

    def CCNOT(self, control1, control2, target):
        return Toffoli(self.psi0, control1, control2, target)

    def mZ(self, target):
        res = Ry(self.psi0, target, np.pi)
        return Rx(self.psi0, target, np.pi)

    def CCCNOT(self, control1, control2, control3, target):
        res = H(self.psi0, target)
        res = CCRz(self.psi0, control2, control3, target, -np.pi/2)
        res = CNOT(self.psi0, control1, control2)
        res = CCRz(self.psi0, control2, control3, target, np.pi/2)
        res = CNOT(self.psi0, control1, control2)
        res = CCRz(self.psi0, control1, control3, target, -np.pi/2)
        res = H(self.psi0, target)
        res = CP(self.psi0, control2, control3, -np.pi/4)
        res = CNOT(self.psi0, control1, control2)
        res = CP(self.psi0, control2, control3, np.pi/4)
        res = CNOT(self.psi0, control1, control2)
        return CP(self.psi0, control1, control3, -np.pi/4)

    def CCCRy(self, control1, control2, control3, target, theta):
        res = CRy(self.psi0, control3, target, theta/2)
        res = CCNOT(self.psi0, control1, control2, control3)
        res = CRy(self.psi0, control3, target, -theta/2)
        res = CCNOT(self.psi0, control1, control2, control3)
        return CCRy(self.psi0, control1, control2, target, theta/2)

    def CCCRz(self, control1, control2, control3, target, theta):
        res = CRz(self.psi0, control3, target, theta/2)
        res = CCNOT(self.psi0, control1, control2, control3)
        res = CRz(self.psi0, control3, target, -theta/2)
        res = CCNOT(self.psi0, control1, control2, control3)
        return CCRz(self.psi0, control1, control2, target, theta/2)

    def CCCP(self, control1, control2, control3, target, theta, b = 0b11):
        res = CP(self.psi0, control3, target, theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
            res = X(self.psi0, control2)
        res = CCNOT(self.psi0, control1, control2, control3)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
            res = X(self.psi0, control2)
        res = CP(self.psi0, control3, target, -theta/2, b = b)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
            res = X(self.psi0, control2)
        res = CCNOT(self.psi0, control1, control2, control3)
        if b == 0b00 or b == 0b01:
            res = X(self.psi0, control1)
            res = X(self.psi0, control2)
        return CCP(self.psi0, control1, control2, target, theta/2, b = b)
