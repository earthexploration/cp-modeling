# -*- coding: utf-8 -*-
import matplotlib
import numpy as np
from scipy.integrate import ode
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import os
sizeOfFont = 16
fontProperties = {'family' : 'serif',
        'serif' : ['Computer Modern Serif'],
        'weight' : 'normal',
        'size' : sizeOfFont
        }
matplotlib.rc('font', **fontProperties)
plt.rc('text', usetex=True)
plt.rc('font', **fontProperties)


# 11, 22, 33, 12, 13, 23

def deps_bulk(dsig, K):
    dJ1 = dsig[0] + dsig[1] + dsig[2]
    
    return dJ1/(3.0*K)

def deps_deve(dsig, mu):
    dJ1 = dsig[0] + dsig[1] + dsig[2]
    dsig[0] -= dJ1/3.0
    dsig[1] -= dJ1/3.0
    dsig[2] -= dJ1/3.0
    
    return dsig/2.0/mu

def f(epsp, sigy, E, N):
    return (1.0 + E/sigy*epsp)**N

def deps_pl(sig, sigy, epsp, deps, E, N, m):
    J1 = sig[0] + sig[1] + sig[2]
    sig[0] -= J1/3.0
    sig[1] -= J1/3.0
    sig[2] -= J1/3.0

    J2 = np.sqrt(3.0/2.0*np.dot(sig, sig))
    
    dj1 = deps[0] + deps[1] + deps[2]
    deps[0] -= dj1/3.0
    deps[1] -= dj1/3.0
    deps[2] -= dj1/3.0

    dj2 = np.sqrt(2.0/3.0*np.dot(deps, deps))

    depsp = dj2*(J2/sigy/f(epsp, sigy, E, N))**m

    return depsp

def deps_devp(sig, sigy, epsp, deps, E, N, m): 
    
    depsp = deps_pl(sig, sigy, epsp, deps, E, N, m)

    J1 = sig[0] + sig[1] + sig[2]
    sig[0] -= J1/3.0
    sig[1] -= J1/3.0
    sig[2] -= J1/3.0

    J2 = np.sqrt(3.0/2.0*np.dot(sig, sig))
    if J2 > 0 :
        deps_p = 3.0*depsp/2.0/J2*sig
        return deps_p
    else:
        return np.zeros(deps.shape)


def g(epsp, y, deps, sigy, E, N, m, nu):
    
    eps_e = y[0:6]
    eps_p = y[6:12]
    #epsp = y[12]

    C = np.zeros((6,6))
    C[0][0] = 1-nu
    C[1][1] = 1-nu
    C[2][2] = 1-nu
    C[0][1] = nu
    C[0][2] = nu
    C[1][2] = nu
    
    C[3][3] = (1-2*nu)/2.0
    C[4][4] = (1-2*nu)/2.0
    C[5][5] = (1-2*nu)/2.0

    for i in range(6):
        for j in range(i):
            C[i][j] = C[j][i]

    C = E/(1+nu)/(1-2*nu)*C

    eps_e[3:6] = 2*eps_e[3:6]
    sig = np.matmul(C, eps_e)
    
    depsp = deps_pl(sig, sigy, epsp, deps, E, N, m)
    deps_p = deps_devp(sig, sigy, epsp, deps, E, N, m)
    deps_e = deps - deps_p
    
    #dydt = y
    #dydt[0:6] = deps_e/depsp
    #dydt[6:12] = deps_p/depsp
    #dydt[12] = depsp
    return np.array([deps_e/depsp,deps_p/depsp]).flatten()

def Hooke(y, E, nu):
    # y = [eps_e, eps_p]
    C = np.zeros((6,6))
    C[0][0] = 1-nu
    C[1][1] = 1-nu
    C[2][2] = 1-nu
    C[0][1] = nu
    C[0][2] = nu
    C[1][2] = nu
    
    C[3][3] = (1-2*nu)/2.0
    C[4][4] = (1-2*nu)/2.0
    C[5][5] = (1-2*nu)/2.0

    for i in range(6):
        for j in range(i):
            C[i][j] = C[j][i]

    C = E/(1+nu)/(1-2*nu)*C
    
    for i in range(len(y)):
        y[i, 6:12] += y[i, 0:6]
        y[i, 9:12] = 2*y[i, 9:12]
        y[i, 3:6] = 2*y[i, 3:6] 
        y[i, 0:6] = np.matmul(C, y[i, 0:6])

    return y

def gfunc(t, y, deps, sigy, E, N, m, nu):
    
    eps_e = y[0:6]
    eps_p = y[6:12]
    epsp = y[12]

    C = np.zeros((6,6))
    C[0][0] = 1-nu
    C[1][1] = 1-nu
    C[2][2] = 1-nu
    C[0][1] = nu
    C[0][2] = nu
    C[1][2] = nu
    
    C[3][3] = (1-2*nu)/2.0
    C[4][4] = (1-2*nu)/2.0
    C[5][5] = (1-2*nu)/2.0

    for i in range(6):
        for j in range(i):
            C[i][j] = C[j][i]

    C = E/(1+nu)/(1-2*nu)*C

    eps_e[3:6] = 2*eps_e[3:6]
    sig = np.matmul(C, eps_e)
    
    depsp = deps_pl(sig, sigy, epsp, deps, E, N, m)
    deps_p = deps_devp(sig, sigy, epsp, deps, E, N, m)
    deps_e = deps - deps_p
    
    dydt = np.zeros(y.shape)
    dydt[0:6] = deps_e    
    dydt[6:12] = deps_p
    dydt[12] = depsp
    return dydt 


def main():

    m = 20
    N = 0.2
    E = 1
    sigy = 0.2E-2*E
    nu = 0.3
    K = E/3.0/(1-2*nu)

    deps = np.array([1e-2, 0, 0, 0, 0, 0])


    t0 = 0
    t1 = 100 #20
    dt = 1e-3
    y0 = np.zeros((13,1))
    r = ode(gfunc,).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0).set_f_params(deps, sigy, E, N, m, nu)
    y = []
    while r.successful() and r.t < t1:
        #print(r.t+dt, r.integrate(r.t+dt))
        y.append(r.integrate(r.t+dt).real)
    y = np.array(y)

   # t2 = 100
   # y0 = y[-1, :]
   # y0[0:6] = 0
   # r.set_initial_value(y0, t1).set_f_params(deps, sigy, E, N, m, nu)
   # y1 = []
   # while r.successful() and r.t < t2:
   #     #print(r.t+dt, r.integrate(r.t+dt))
   #     y1.append(r.integrate(r.t+dt).real)
   # y1 = np.array(y1)

    C = np.zeros((6,6))
    C[0][0] = 1-nu
    C[1][1] = 1-nu
    C[2][2] = 1-nu
    C[0][1] = nu
    C[0][2] = nu
    C[1][2] = nu
    
    C[3][3] = (1-2*nu)/2.0
    C[4][4] = (1-2*nu)/2.0
    C[5][5] = (1-2*nu)/2.0

    for i in range(6):
        for j in range(i):
            C[i][j] = C[j][i]

    C = E/(1+nu)/(1-2*nu)*C

    for i in range(len(y)):
        y[i, 6:12] += y[i, 0:6]
        y[i, 9:12] = 2*y[i, 9:12]
        y[i, 3:6] = 2*y[i, 3:6] 
        y[i, 0:6] = np.matmul(C, y[i, 0:6])
    
    #Hooke(y1, E, nu)
    
    #sig0 = sigy*f(epsp, sigy, E, N)
    #plt.plot(epsp+sig0/E, sig0/sigy)

    plt.plot(y[:, 6], y[:, 0]/sigy, linestyle='-')
    #plt.plot(y1[:, 6], y1[:, 0]/sigy, linestyle='--')
    plt.xlabel(r'$\varepsilon_{11}$')
    plt.ylabel(r'$\sigma/\sigma_{y}$')
    plt.grid() 
    plt.tight_layout()
    plt.show()
    return 0

main()

