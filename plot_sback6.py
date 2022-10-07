import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.integrate import odeint
import os

def dot_voce(s, p, dot_p, M, h0, tau0, taus):

    c0 = tau0
    c1 = (taus-tau0)
    c2 = h0/(taus-tau0) 
    # 
    dot_s = (c0 + c1)*c2*abs(dot_p) - c2*abs(dot_p)*s
    
    return dot_s

def dot_sback(sb, p, dot_p, M, h0, tau0, taus, hs, g):
    # p: equivalent plastic strain

    F = (h0-hs)*math.pow(1.0/math.cosh((h0-hs)*p/(taus-tau0)), 2.0)+hs
    g0 = (taus-tau0)/h0
    G = 1.0 + M*math.tanh(g/g0)
    h = F*G

    dot_sb = h*dot_p -g*abs(dot_p)*sb

    return dot_sb

def f(y, t, dot_p, M, h0, tau0, taus, hs, g):

    p = y[0]
    sb = y[1]
    sv = y[2]

    dot_sb = dot_sback(sb, p, dot_p, M, h0, tau0, taus, hs, g)
    dot_sv = dot_voce(sv, p, dot_p, M, h0, tau0, taus)
   
    #dydt = [abs(dot_p), dot_sb, dot_sv]
    dydt = [dot_p, dot_sb, dot_sv]

    return dydt

def main():
    M = 2.0 #4
    h0 = 56300 #610 #5000 #150000
    tau0 = 702 #100 #50
    taus = 1510 #235 #100 #200
    g = (taus-tau0) #3000
    hs = 1.4*h0 #150000 #150000 #50000


    dot_p = 1e-3
    #y0 = [0, 0, M*np.sign(dot_p)*tau0]
    y0 = [0, 0, np.sign(dot_p)*tau0]
    #y0 = [0, 0, 0]

    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g))
    
    p1 = sol[:,0]
    s1 = sol[:,1] + np.sign(dot_p)*sol[:,2]

    plt.plot(p1, s1, '-', label="prestrain")

    #print(sol[:,1])
    dot_p = 1e-3
    y0 = sol[-1,:]
    t = np.linspace(0, 150, 501)
    sol0 = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g)) 
    p0 = sol0[:,0]
    s0 = sol0[:,1] + np.sign(dot_p)*sol0[:,2]
    #plt.plot(p0, s0, '-')
    plt.plot(y0[0]+abs(p0-y0[0]), abs(s0),'--', label="monotonic")


    dot_p = -1e-3
    y0 = sol[-1,:]
    #y0[2] = M*np.sign(dot_p)*tau0 
    #y0[2] = sol[-1,2]
    #y0[1] = np.sign(dot_p)*sol[-1,1]
    #y0[2] = np.sign(dot_p)*sol[-1,2]
    #y0[2] = np.sign(dot_p)*s1[-1]
    y0[2] = abs(s1[-1])

       
    
    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g))
    
    p2 = sol[:,0]
    s2 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    #plt.plot(y0[0]+abs(p2-y0[0]), np.sign(dot_p)*s2)
    plt.plot(y0[0]+abs(p2-y0[0]), abs(s2),':', label="1st compression")
    #plt.plot(p2, s2, '--')
    #y1 = y0[0]+abs(p2[-1]-y0[0])
    
    dot_p = 1e-3
    y0 = sol[-1,:]
    #y0[1] = np.sign(dot_p)*sol[-1,1]
    #y0[2] = sol[-1,2] 
    #y0[2] = np.sign(dot_p)*s2[-1]
    y0[2] = abs(s2[-1])
     
    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g))
    
    p3 = sol[:,0]
    s3 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    #plt.plot(p3, s3, '-')
    plt.plot(y0[0]+abs(p3-y0[0]), abs(s3),'-.', label="1st tension")
    
    dot_p = -1e-3
    y0 = sol[-1,:]
    #y0[1] = np.sign(dot_p)*sol[-1,1] 
    #y0[2] = sol[-1,2] 
    #y0[2] = np.sign(dot_p)*s3[-1]
    y0[2] = abs(s3[-1])

    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g))
    
    p4 = sol[:,0]
    s4 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    #plt.plot(p4, s4, '--')
    plt.plot(y0[0]+abs(p4-y0[0]), abs(s4),':', label="2nd compression")
    #y1 = y0[0]+abs(p4[-1]-y0[0])
    
    dot_p = 1e-3
    y0 = sol[-1,:]
    y0[2] = abs(s4[-1])
    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g)) 
    p5 = sol[:,0]
    s5 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    #plt.plot(p5, s5, '-')
    plt.plot(y0[0]+abs(p5-y0[0]), abs(s5),'-.', label="2nd tension")
    
    dot_p = -1e-3
    y0 = sol[-1,:]
    y0[2] = abs(s5[-1])
    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g)) 
    p6 = sol[:,0]
    s6 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    #plt.plot(p6, s6, '--')
    plt.plot(y0[0]+abs(p6-y0[0]), abs(s6),':', label="3rd compression") 
    #y1 = y0[0]+abs(p6[-1]-y0[0])

    dot_p = 1e-3
    y0 = sol[-1,:]
    y0[2] = abs(s6[-1])
    t = np.linspace(0, 100, 501)
    sol = odeint(f, y0, t, args=(dot_p, M, h0, tau0, taus, hs, g)) 
    p7 = sol[:,0]
    s7 = sol[:,1] + np.sign(dot_p)*sol[:,2]
    plt.plot(p7, s7, '-') 
    plt.plot(y0[0]+abs(p7-y0[0]), abs(s7),'-.', label="3rd tension")

    plt.xlabel(r"$\varepsilon_{p}$")
    plt.ylabel(r"$\sigma$")
    plt.title(r"Voce's law: $\sigma(\varepsilon_{p})$ ")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


main()
