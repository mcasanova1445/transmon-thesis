# Importando librerias y modulos necesarios

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from qutip import *

def gaussianpulse(x,ts,tf):
    s = (tf-ts)/6
    m = (ts+tf)/2
    return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))*norm.pdf(x, loc = m, scale = s)/0.997300204

def squarepulse(x,ts,tf):
    s = (tf-ts)/6
    m = (ts+tf)/2
    return (np.heaviside(x-m+3*s,1)-np.heaviside(x-m-3*s,1))/(6*s)

def plot_wigner_2d(res):
    #fig, axes = plt.subplots(1, 2, subplot_kw={'projection': '3d'}, figsize=(12, 6))
    fig = plt.figure(figsize=(8, 8))

    ax = fig.add_subplot(1, 1, 1)
    plot_wigner(ptrace(res.states[-1],0), fig=fig, ax=ax, alpha_max=10);

    fig.tight_layout()

def plot_drive_expect(res,args):
    tlist = res.times

    if args == 0:
        fig, axes = plt.subplots(2, 1, sharex=True, figsize=(12,8))

        axes[0].plot(tlist, np.real(expect(n, res.states)), 'r', linewidth=2, label="cavity")
        axes[0].set_ylabel("Occupation probability", fontsize=16)
        axes[0].legend()

        axes[1].plot(tlist, np.real(expect(qop('n',0), res.states)), 'b', linewidth=2, label="qubit 0")
        axes[1].plot(tlist, np.real(expect(qop('n',1), res.states)), 'g', linewidth=2, label="qubit 1")
        axes[1].plot(tlist, np.real(expect(qop('n',2), res.states)), 'c', linewidth=2, label="qubit 2") #MELCSCELDQ
        axes[1].plot(tlist, np.real(expect(qop('n',3), res.states)), 'm', linewidth=2, label="qubit 3") #MELCSCELDQ
        # axes[1].plot(tlist, np.real(expect(qop('n',4), res.states)), 'r', linewidth=2, label="qubit 4") #MELCSCELDQ
        # axes[1].plot(tlist, np.real(expect(qop('n',5), res.states)), 'y', linewidth=2, label="qubit 5") #MELCSCELDQ
        # axes[1].plot(tlist, np.real(expect(qop('n',6), res.states)), 'k', linewidth=2, label="qubit 6") #MELCSCELDQ
        # axes[1].plot(tlist, np.real(expect(qop('n',7), res.states)), 'w', linewidth=2, label="qubit 7") #MELCSCELDQ
        axes[1].set_ylim(0, 1)

        axes[1].set_xlabel("Time (ns)", fontsize=16)
        axes[1].set_ylabel("Occupation probability", fontsize=16)
        axes[1].legend()

    else:
        fig, axes = plt.subplots(3, 1, sharex=True, figsize=(12,12))

        axes[0].plot(tlist, np.array(list(ksi_t(tlist,args))) / (2*np.pi), 'b', linewidth=2, label="drive envelope")
        #axes[0].set_ylim(-0.3, 0.3) # *np.cos(wd*(tlist))
        axes[0].set_ylabel("Energy (GHz)", fontsize=16)
        axes[0].legend()

        axes[1].plot(tlist, np.real(expect(n, res.states)), 'r', linewidth=2, label="cavity")
        axes[1].set_ylabel("Occupation probability", fontsize=16)
        axes[1].legend()

        axes[2].plot(tlist, np.real(expect(qop('n',0), res.states)), 'b', linewidth=2, label="qubit 0")
        axes[2].plot(tlist, np.real(expect(qop('n',1), res.states)), 'g', linewidth=2, label="qubit 1")
        axes[2].plot(tlist, np.real(expect(qop('n',2), res.states)), 'c', linewidth=2, label="qubit 2") #MELCSCELDQ
        axes[2].plot(tlist, np.real(expect(qop('n',3), res.states)), 'm', linewidth=2, label="qubit 3") #MELCSCELDQ
        # axes[2].plot(tlist, np.real(expect(qop('n',4), res.states)), 'r', linewidth=2, label="qubit 4") #MELCSCELDQ
        # axes[2].plot(tlist, np.real(expect(qop('n',5), res.states)), 'y', linewidth=2, label="qubit 5") #MELCSCELDQ
        # axes[2].plot(tlist, np.real(expect(qop('n',6), res.states)), 'k', linewidth=2, label="qubit 6") #MELCSCELDQ
        # axes[2].plot(tlist, np.real(expect(qop('n',7), res.states)), 'w', linewidth=2, label="qubit 7") #MELCSCELDQ
        axes[2].set_ylim(0, 1)

        axes[2].set_xlabel("Time (ns)", fontsize=16)
        axes[2].set_ylabel("Occupation probability", fontsize=16)
        axes[2].legend()

    fig.tight_layout()

# Parametros del sistema

N = 50

wr = 10.0 * 2 * np.pi
# wq = np.array([5.0 * 2 * np.pi, 6.0 * 2 * np.pi, 7.0 * 2 * np.pi, 8.0 * 2 * np.pi])
wq = np.array([5.0 * 2 * np.pi, 6.0 * 2 * np.pi, 7.0 * 2 * np.pi, 8.0 * 2 * np.pi, 11.0 * 2 * np.pi, 12.0 * 2 * np.pi, 13.0 * 2 * np.pi, 14.0 * 2 * np.pi])
wq_swap = 9 * 2 * np.pi

g = np.array([0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi, 0.1 * 2*np.pi])

D = wq - wr

chi = g**2 / abs(wr-wq)

kappa = 0.0001
gamma = np.array([0.005, 0.005, 0.005, 0.005])

'''
for i in range(4):
    for j in range(4):
        L[i,j] = g[i]*g[j]*((1/D[i]) + (1/D[j]))/2

L = np.matrix(L)
'''

# cavity operators
a = tensor(destroy(N), qeye(2), qeye(2), qeye(2), qeye(2)) # Modificar esta linea cuando se cambie el numero de qubits #MELCSCELDQ
# a = tensor(destroy(N), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2), qeye(2)) # Modificar esta linea cuando se cambie el numero de qubits #MELCSCELDQ
n = a.dag() * a
Id = tensor(qeye(N), qeye(2), qeye(2), qeye(2), qeye(2))

def qop_part(operator, target):
    if target == 0:
        qop_dict = {'sm' : destroy(2), 'sp' : (destroy(2)).dag(), 
                    'sx' : sigmax(), 'sy' : sigmay(), 'sz' : sigmaz(),
                    'n' : (destroy(2)).dag() * destroy(2)}
        return qop_dict[operator]
    else:
        return qeye(2)

def qop(operator, target):
    return tensor(qeye(N), qop_part(operator, target-0), qop_part(operator, target-1), qop_part(operator, target-2), qop_part(operator, target-3)) #MELCSCELDQ
    
# def qop(operator, target):
    # return tensor(qeye(N), qop_part(operator, target-0), qop_part(operator, target-1), qop_part(operator, target-2), qop_part(operator, target-3), qop_part(operator, target-4), qop_part(operator, target-5), qop_part(operator, target-6), qop_part(operator, target-7)) #MELCSCELDQ

c_ops = [np.sqrt(kappa) * a, np.sqrt(gamma[0]) * qop('sm', 0), np.sqrt(gamma[1]) * qop('sm', 1), np.sqrt(gamma[2]) * qop('sm', 2), np.sqrt(gamma[3]) * qop('sm', 3)]

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

def Rx(psi0, target, theta):
    tlist = np.linspace(0, 10, 5000)

    wd = wq[target]

    Dr = wr-wd
    Dq = wq-wd

    Hsyst = Dr*n #+ Dq[0]*qop('sz',0)/2 + Dq[1]*qop('sz',1)/2 + Dq[2]*qop('sz',2)/2 + Dq[3]*qop('sz',3)/2 #MELCSCELDQ
    for i in range(4): #MELCSCELDQ
        Hsyst = Hsyst - Dq[i]*qop('sz',i)/2

    # Hsyst = Dr*n + Dq[0]*qop('sz',0)/2 + Dq[1]*qop('sz',1)/2 + Dq[2]*qop('sz',2)/2 + Dq[3]*qop('sz',3)/2 + Dq[4]*qop('sz',4)/2 + Dq[5]*qop('sz',5)/2 + Dq[6]*qop('sz',6)/2 + Dq[7]*qop('sz',7)/2 #MELCSCELDQ
    H_t = [[qop('sx',target)/2, ksi_t], Hsyst]

    args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : wq[target]}
    res = mesolve(H_t, psi0, tlist, [], [], args = args)

    # plot_drive_expect(res,args)

    return res

def Ry(psi0, target, theta):
    tlist = np.linspace(0, 10, 5000)

    wd = wq[target]

    Dr = wr-wd
    Dq = wq-wd

    Hsyst = Dr*n #+ Dq[0]*qop('sz',0)/2 + Dq[1]*qop('sz',1)/2 + Dq[2]*qop('sz',2)/2 + Dq[3]*qop('sz',3)/2 #MELCSCELDQ
    for i in range(4): #MELCSCELDQ
        Hsyst = Hsyst - Dq[i]*qop('sz',i)/2

    H_t = [[qop('sy',target)/2, ksi_t], Hsyst]

    args = {'A' : theta, 'ts' : 0, 'tf' : 10, 'w' : wq[target]}
    res = mesolve(H_t, psi0, tlist, [], [], args = args)

    # plot_drive_expect(res,args)

    return res

def Rz(psi0, target, theta):
    res = Ry(psi0, target, np.pi/2)
    res = Rx(res.states[-1], target, theta)
    return Ry(res.states[-1], target, -np.pi/2)

def X(psi0, target):
    return Rx(psi0, target, np.pi)

def Y(psi0, target):
    return Ry(psi0, target, np.pi)

def Z(psi0, target, theta):
    res = Rx(psi0, target, np.pi/2)
    res = Y(res.states[-1], target)
    return Rx(res.states[-1], target, -np.pi/2)

def H(psi0, target):
    res = Ry(psi0, target, np.pi/2)
    return X(res.states[-1], target)

def sqrtiSWAP(psi0, target1, target2):
    wqt1 = wq[target1]
    wq[target1] = wq_swap
    
    wqt2 = wq[target2]
    wq[target2] = wq_swap

    D = wq - wr

    # J = np.abs(g[target1] * g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2])) # Revisar Hamiltoniano #RH
    J = np.abs(g[target1] * g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

    tf = np.pi/(4*J)
    tlist = np.linspace(0, tf, 5000)

    Hsyst = wr*n #+ Dq[0]*qop('sz',0)/2 + Dq[1]*qop('sz',1)/2 + Dq[2]*qop('sz',2)/2 + Dq[3]*qop('sz',3)/2 #MELCSCELDQ
    for i in range(4): #MELCSCELDQ
        # Hsyst = Hsyst - wq[i]*(Id/2 + 2 * (g[i]**2/D[i])*(n + Id/2))*qop('sz',i) #RH
        Hsyst = Hsyst - wq[i]*qop('sz',i)/2 #RH

    for i in range(4):
        for j in range(4):
            if i!=j:
                # Hsyst = Hsyst + (g[i]*g[j]/D[j])*(qop('sp',i)*qop('sm',j) + qop('sm',i)*qop('sp',j)) #RH
                Hsyst = Hsyst + (g[i]*g[j]/D[j])*(qop('sp',i)*qop('sm',j) + qop('sm',i)*qop('sp',j))/2 #RH

    res = mesolve(Hsyst, psi0, tlist, [], [])

    wq[target1] = wqt1
    wq[target2] = wqt2
    D = wq - wr

    args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]}
    # plot_drive_expect(res,args)

    return res

def iSWAP(psi0, target1, target2):
    wqt1 = wq[target1]
    wq[target1] = wq_swap
    
    wqt2 = wq[target2]
    wq[target2] = wq_swap

    D = wq - wr

    # J = np.abs(g[target1] * g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2])) # Revisar Hamiltoniano #RH
    J = np.abs(g[target1] * g[target2] * (D[target1] + D[target2]) / (D[target1] * D[target2]))/2 # Revisar Hamiltoniano #RH

    tf = 0.9*np.pi/(2*J)
    tlist = np.linspace(0, tf, 2*5000)

    Hsyst = wr*n #+ Dq[0]*qop('sz',0)/2 + Dq[1]*qop('sz',1)/2 + Dq[2]*qop('sz',2)/2 + Dq[3]*qop('sz',3)/2 #MELCSCELDQ
    for i in range(4): #MELCSCELDQ
        # Hsyst = Hsyst - wq[i]*(Id/2 + 2 * (g[i]**2/D[i])*(n + Id/2))*qop('sz',i) #RH
        Hsyst = Hsyst - wq[i]*qop('sz',i)/2 #RH

    for i in range(4):
        for j in range(4):
            if i!=j:
                # Hsyst = Hsyst + (g[i]*g[j]/D[j])*(qop('sp',i)*qop('sm',j) + qop('sm',i)*qop('sp',j)) #RH
                Hsyst = Hsyst + (g[i]*g[j]/D[j])*(qop('sp',i)*qop('sm',j) + qop('sm',i)*qop('sp',j))/2 #RH

    res = mcsolve(Hsyst, psi0, tlist, c_ops, [])

    wq[target1] = wqt1
    wq[target2] = wqt2
    D = wq - wr

    args = {'A' : 0, 'ts' : 0, 'tf' : tf, 'w' : wq[target1]}
    plot_drive_expect(res,args)

    return res

def CNOT(psi0, control, target):
    res = H(psi0, target)
    res = Rz(res.states[-1], target, -np.pi/2)
    res = Rz(res.states[-1], control, -np.pi/2)
    res = iSWAP(res.states[-1], control, target)
    res = H(res.states[-1], control)
    res = iSWAP(res.states[-1], control, target)
    res = Rx(res.states[-1], target, np.pi/2)
    res = iSWAP(res.states[-1], control, target)
    res = Rx(res.states[-1], control, np.pi/2)
    res = iSWAP(res.states[-1], control, target)
    return Rx(res.states[-1], target, np.pi/2)

def CRy(psi0, control, target, theta): #Las rotaciones aqui son al reves que las que tengo en Mathematica
    res = Ry(psi0,target,theta/2)
    res = CNOT(res.states[-1],control,target)
    res = Ry(res.states[-1],target,-theta/2)
    return CNOT(res.states[-1],control,target)

def CRz(psi0, control, target, theta): #Las rotaciones aqui son al reves que las que tengo en Mathematica
    res = Rz(psi0,target,theta/2)
    res = CNOT(res.states[-1],control,target)
    res = Rz(res.states[-1],target,-theta/2)
    return CNOT(res.states[-1],control,target)

def Dis(psi0):
    tf = 0.05
    tlist = np.linspace(0, tf, 5000)

    wd = wr

    Hsyst = wr*(a.dag()*a + Id/2)
    for i in range(4): #MELCSCELDQ
        Hsyst = Hsyst - wq[i]*qop('sz',i)/2 + g[i]*qop('sx',i)*(a+a.dag())

    
    H_t = [[a, ksi_tm],[a.dag(),ksi_tp], Hsyst]

    args = {'A' : wr/(3*np.pi), 'ts' : 0, 'tf' : tf, 'w' : wd}
    res = mesolve(H_t, psi0, tlist, [], [], args = args)

    # plot_drive_expect(res,args)

    plot_wigner_2d(res)

    return res

def Dis_wait(psi0):
    tlist = np.linspace(0, 2*np.pi/chi[0], 50000)

    Hsyst = wr*(a.dag()*a + Id/2)
    for i in range(4): #MELCSCELDQ
        Hsyst = Hsyst - wq[i]*qop('sz',i)/2 + g[i]*qop('sx',i)*(a+a.dag())

    for i in range(4):
        for j in range(4):
            Hsyst = Hsyst + (g[i]*g[j]*((1/D[i]) + (1/D[j]))/2)*qop('sx',i)*qop('sx',j)/2

    res = mesolve(Hsyst, psi0, tlist, [], [])

    # plot_drive_expect(res,0)

    plot_wigner_2d(res)

    return res
