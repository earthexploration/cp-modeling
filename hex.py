#!/usr/python
import math
import numpy as np
import scipy
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from matplotlib import cm
import matplotlib.colors as colors
import matplotlib
matplotlib.use('Qt5Agg')
#from matplotlib.animation import ImageMagickWriter
sizeOfFont = 12
fontProperties = {'family' : 'serif',
        'serif' : ['COmputer Modern Serif'],
        'weight' : 'normal',
        'size' : sizeOfFont
        }
matplotlib.rc('font', **fontProperties)
plt.rc('text', usetex=True)
plt.rc('font', **fontProperties)
class Hex:
  a1 = [-0.5, -0.5*math.sqrt(3.0), 0.0,]
  a2 = [1.0, 0.0, 0.0,]
  a3 = [-0.5, 0.5*math.sqrt(3.0), 0.0,]
  c = [0.0, 0.0, math.sqrt(8.0/3.0),]

  g = [a1, a2, a3, c]

  nslip = [3, 3, 3, 6, 12, 6]
  slipdir13 = [[ 2, -1, -1, 0,],
               [-1,  2, -1, 0,],
               [-1, -1,  2, 0,],]

  slipnor13 = [[0, 0, 0, 1,],
               [0, 0, 0, 1,],
               [0, 0, 0, 1,],]


  slipdir23 = [[ 2, -1, -1, 0,],
               [-1,  2, -1, 0,],
               [-1, -1,  2, 0,],]

  slipnor23 = [[ 0,  1,  -1,  0,],
               [-1,  0,   1,  0,],
               [1,  -1,  0,  0,],]


  slipdir33 = [[ 0,  1, -1, 0,],
               [-1,  0,  1, 0,],
               [ 1, -1,  0, 0,],]

  slipnor33 = [[ 2, -1,  -1, 0,],
               [-1,  2,  -1, 0,],
               [-1, -1,   2, 0,],]


  slipdir46 = [[ 2, -1, -1, 0,],
               [-1, 2, -1, 0,],
               [-1, -1, 2, 0,],
               [1, 1, -2, 0,],
               [-2, 1, 1, 0,],
               [1, -2, 1, 0,],]

  slipnor46 = [[0, 1, -1, 1,],
               [-1, 0, 1, 1,],
               [1, -1, 0, 1,],
               [-1, 1, 0, 1,],
               [0, -1, 1, 1,],
               [1, 0, -1, 1,],]

  slipdir512 = [[2, -1, -1, 3,],
                [1, -2, 1, 3,],
                [-1, -1, 2, 3,],
                [-2, 1, 1, 3,],
                [-1, 2, -1, 3,],
                [1, 1, -2, 3,],
                [-2, 1, 1, 3,],
                [-1, 2, -1, 3,],
                [1, 1, -2, 3,],
                [2, -1, -1, 3,],
                [1, -2, 1, 3,],
                [-1, -1, 2, 3,],]

  slipnor512 = [[ -1, 1, 0, 1,],
                [-1, 1, 0, 1,],
                [1, 0, -1, 1],
                [1, 0, -1, 1,],
                [0, -1, 1, 1,],
                [0, -1, 1, 1,],
                [1, -1, 0, 1,],
                [1, -1, 0, 1,],
                [-1, 0, 1, 1,],
                [-1, 0, 1, 1,],
                [0, 1, -1, 1,],
                [0, 1, -1, 1,],]


  slipdir66 = [[2, -1, -1, 3,],
               [-1, 2, -1, 3,],
               [-1, -1, 2, 3,],
               [-2, 1, 1, 3,],
               [1, -2, 1, 3,],
               [1, 1, -2, 3,],]

  slipnor66 = [[-2, 1, 1, 2,],
               [1, -2, 1, 2,],
               [1, 1, -2, 2,],
               [2, -1, -1, 2,],
               [-1, 2, -1, 2,],
               [-1, -1, 2, 2,],]

  fcc_slipnor = [[1,  1, -1,], # primary system
                 [1, -1,  1,], # conjugate system
                 [1, -1, -1,], # cross-glide system
                 [1,  1,  1,],] # critical system

  # primary system
  fcc_slipnor_primary = [1, 1, -1,]
  fcc_slipdir_primary = [[ 1,  0,  1,],
                         [ 0, -1, -1,],
                         [-1,  1,  0,],]
  # conjuate system
  fcc_slirnor_conjugate = [1, -1, 1,]
  fcc_slipdir_conjugate = [[-1, -1,  0],
                           [ 0,  1,  1],
                           [ 1,  0, -1],]

  # cross-glide system
  fcc_slipnor_crossglide = [1, -1, -1,]
  fcc_slipdir_crossglide = [[-1, 0, -1],
                            [ 1, 1,  0],
                            [ 0, -1, 1],]

  # critical system
  fcc_slipnor_critical = [1, 1, 1,]
  fcc_slipdir_critical = [[-1, 0,  1],
                          [ 0, 1, -1],
                          [ 1, -1, 0],]




  c11 = 160.0e6
  c12 = 90.0e6
  c13 = 66.0e6
  c33 = 181.7e6
  c44 = 46.5e6

  C = [[c11, c12, c13, 0, 0, 0,],
       [c12, c11, c13, 0, 0, 0,],
       [c13, c13, c33, 0, 0, 0,],
       [0, 0, 0, c44, 0, 0,],
       [0, 0, 0, 0, c44, 0,],
       [0, 0, 0, 0, 0, (c11-c12)/2.0]]

  def rot_ZXZ(alpha, beta, gamma):
    R = [[math.cos(alpha)*math.cos(gamma)-math.cos(beta)*math.sin(alpha)*math.sin(gamma),
          math.sin(alpha)*math.cos(gamma)+math.cos(beta)*math.cos(alpha)*math.sin(gamma),
          math.sin(beta)*math.sin(gamma),],
         [-math.cos(alpha)*math.sin(gamma)-math.cos(beta)*math.sin(alpha)*math.cos(gamma),
          -math.sin(alpha)*math.sin(gamma)+math.cos(beta)*math.cos(alpha)*math.cos(gamma),
          math.sin(beta)*math.cos(gamma),],
         [math.sin(beta)*math.sin(alpha), -math.sin(beta)*math.cos(alpha), math.cos(beta),],]
    return R

  n = 10
  a = 0.001

  def hard_powerlaw_strainrate(gamma, tauslp, gslip, slipnor, slipdir, prop):
    n = prop[0]
    a = prop[1]
    x = tauslp/gslip

    fslip = a*abs(x)**(n)*np.sign(x)
    return fslip

  def hard_powerlaw_dfdxsp(gamma, tauslp, gslip, slipnor, slipdir, prop):
    n = prop[0]
    a = prop[1]
    x = tauslp/gslip

    dfdxsp = a*n*abs(x)**(n-1)
    return dfdxsp

  h0 = 541.5
  tau_s = 109.5
  tau_0 = 60.8

  def hard_moduli_asaro1983(gamma, gmsltl, gamtol, slipnor, slipdir, nslptl, prop):
    h0 = prop[0]
    tau_s = prop[1]
    tau_0 = prop[2]

    term1 = abs(h0*gamtol/(tau_s-tau_0))
    term2 = 1.0/math.cosh(term1)
    hself = h0*math.pow(1.0/math.cosh(term1),2.0)

    q = np.ones(len(gmsltl))*math.sqrt(2)
    q[nslptl] = 1.0
    hlatnt = hself*q
    return hlatnt

  def hard_moduli_asaro1983_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nslptl, prop):
    h0 = prop[0]
    tau_s = prop[1]
    tau_0 = prop[2]

    term1 = abs(h0*gamtol/(tau_s-tau_0))
    term2 = 1.0/math.cosh(term1)
    term3 = h0/(tau_s-tau_0)*np.sign(gamma)
    dhself = -2.0*h0*math.pow(term2,2.0)*math.tanh(term1)*term3

    q = np.ones(len(gmsltl))*math.sqrt(2)
    q[nslptl] = 1.0
    term3 = h0/(tau_s-tau_0)*np.sign(gmsltl)*q
    dhlatnt = -2.0*h0*math.pow(term2,2.0)*math.tanh(term1)*term3
    return dhlatnt

  hs = 15
  gamma0 = 0.001
  def fab(nslipnor1, nslipdir1, nslipnor2, nslipdir2, a1, a2, a3, a4, a5):
    slipdir1 = [0, 0, 0]
    slipdir2 = [0, 0, 0]

    if nslipnor1 == 0:
      slipdir1 = Hex.fcc_slipdir_primary[nslipdir1]
    if nslipnor1 == 1:
      slipdir1 = Hex.fcc_slipdir_conjugate[nslipdir1]
    if nslipnor1 == 2:
      slipdir1 = Hex.fcc_slipdir_crossglide[nslipdir1]
    if nslipnor1 == 3:
      slipdir1 = Hex.fcc_slipdir_critical[nslipdir1]

    if nslipnor2 == 0:
      slipdir2 = Hex.fcc_slipdir_primary[nslipdir2]
    if nslipnor2 == 1:
      slipdir2 = Hex.fcc_slipdir_conjugate[nslipdir2]
    if nslipnor2 == 2:
      slipdir2 = Hex.fcc_slipdir_crossglide[nslipdir2]
    if nslipnor2 == 3:
      slipdir2 = Hex.fcc_slipdir_critical[nslipdir2]

    if nslipnor1 != nslipnor2:
      if np.dot(slipdir1, slipdir2) == -2:
        #print('N: no junction')
        return a1
      if np.dot(slipdir1, slipdir2) == 0:
        #print('H: Hirth lock')
        return a2
      if np.dot(slipdir1, slipdir2) == 1:
        #print('G: glissile junction')
        return a4
      if np.dot(slipdir1, slipdir2) == -1:
        #print('S: sessile junction')
        return a5

    if nslipnor1 == nslipnor2:
      if np.dot(slipdir1, slipdir2) == -1:
        #print('C: coplanar junction')
        return a3

    return 0

  def fcc_fab(a1, a2, a3, a4, a5):
    fcc_fab = np.ones((4, 3, 4, 3))
    for i in range(4):
      for j in range(3):
        for k in range(4):
          for l in range(3):
            #fcc_fab[i, j, k, l] = Hex.fab(i, j, k, l, 8, 8, 8, 15, 20)
            #fcc_fab[i, j, k, l] = Hex.fab(i, j, k, l, 1, 2, 3, 4, 5)
            fcc_fab[i, j, k, l] = Hex.fab(i, j, k, l, a1, a2, a3, a4, a5)
    return fcc_fab.reshape(12, 12)

  def hard_moduli_Bassani1991(gamma, gmsltl, gamtol, slipnor, slipdir, nslptl, prop):
    h0 = prop[0]
    tau_I = prop[1]
    tau_0 = prop[2]
    hs = prop[3]
    gamma0 = prop[4]

    a1 = a2 = a3 = prop[5]
    a4 = prop[6]
    a5 = prop[7]

    #fab = Hex.fcc_fab(a1, a2, a3, a4, a5)[3*nslipnor+nslipdir]
    fab = Hex.fcc_fab(a1, a2, a3, a4, a5)[nslptl]

    term1 = (h0-hs)*gamma/(tau_I-tau_0)
    term2 = 1.0/math.cosh(term1)
    F = (h0-hs)*math.pow(term2,2.0)+hs
    G = 1
    for i, gmsl in enumerate(gmsltl):
      if i != nslptl:
        G = G + fab[i]*math.tanh(gmsl/gamma0)
    hself = F*G
    Q = np.ones(len(gmsltl))*math.sqrt(2)
    Q[nslptl] = 1.0
    hlatn = hself*Q

    return hlatn

  def hard_moduli_Bassani1991_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nslptl, prop):
    h0 = prop[0]
    tau_I = prop[1]
    tau_0 = prop[2]
    hs = prop[3]
    gamma0 = prop[4]

    a1 = a2 = a3 = prop[5]
    a4 = prop[6]
    a5 = prop[7]

    #fab = Hex.fcc_fab(a1, a2, a3, a4, a5)[3*nslipnor+nslipdir]
    fab = Hex.fcc_fab(a1, a2, a3, a4, a5)[nslptl]

    term1 = (h0-hs)*gamma/(tau_I-tau_0)
    term2 = 1.0/math.cosh(term1)
    term3 = (h0-hs)/(tau_I-tau_0)
    #F = (h0-hs)*math.pow(term2,2.0)+hs
    F = -2.0*(h0-hs)*math.pow(term2,2.0)*math.tanh(term1)*term3
    G = 1
    for i, gmsl in enumerate(gmsltl):
      if i != nslptl:
        G = G + fab[i]*math.tanh(gmsl/gamma0)
    dhself = np.zeros(len(gmsltl))
    dhlatn = np.zeros((len(gmsltl),len(gmsltl)))
    dhself[nslptl] = F*G
    dhlatn[nslptl, nslptl] = F*G

    F = (h0-hs)*math.pow(term2,2.0)+hs
    for i, gmsl in enumerate(gmsltl):
      if i != nslptl:
        term4 = gmsl/gamma0
        term5 = 1.0/math.cosh(term4)
        G = fab[i]/gamma0*math.pow(term5,2.0)
        dhself[i] = F*G
        dhlatn[nslptl, i] = F*G
    Q = np.ones(len(gmsltl))*math.sqrt(2)
    Q[nslptl] = 1.0
    for j, q in enumerate(Q):
      if j != nslptl:
        dhlatn[j] = q*dhlatn[nslptl]

    return dhlatn







#R = Hex.rot_ZXZ(257.5, 53.3, 157.3)
#print(np.linalg.norm(Hex.g, axis=1))
csys = np.transpose(Hex.g)
slipdir = np.matmul(np.transpose(Hex.g), np.transpose(Hex.slipdir512))
slipdir = slipdir/np.linalg.norm(slipdir, axis=0)
slipnor = np.matmul(np.transpose(Hex.g), np.transpose(Hex.slipnor512))
slipnor = slipnor/np.linalg.norm(slipnor, axis=0)
#print(slipnor)
#print(slipdir)
#fcc_fab = np.ones((4, 3, 4, 3))
#for i in range(4):
#  for j in range(3):
#    for k in range(4):
#      for l in range(3):
#        #fcc_fab[i, j, k, l] = Hex.fab(i, j, k, l, 8, 8, 8, 15, 20) 
#        fcc_fab[i, j, k, l] = Hex.fab(i, j, k, l, 1, 2, 3, 4, 5) 
#print(fcc_fab.reshape(12, 12))
print(Hex.fcc_fab(1, 2, 3, 4, 5))
nstart = np.sum(Hex.nslip[0:5])
nstop = nstart+Hex.nslip[4]
gmsltl = Hex.a*np.zeros(Hex.nslip[4])
gamtol = np.sum(np.abs(gmsltl))
prop = [Hex.h0, Hex.tau_s, Hex.tau_0]
t = np.arange(0, 0.5, Hex.a)
hself = np.zeros(len(t))
dhself = np.zeros(len(t))
hlatnt = np.zeros((len(gmsltl), len(t)))
dhlatnt = np.zeros((len(gmsltl), len(t)))
nactive = 0

prop1 = [Hex.h0, Hex.tau_s, Hex.tau_0, Hex.hs, Hex.gamma0, 8, 15, 20]
hself1 = np.zeros(len(t))
dhself1 = np.zeros((len(gmsltl), len(t)))
hlatn = np.zeros((len(gmsltl), len(t)))
dhlatn = np.zeros((len(gmsltl), len(gmsltl), len(t)))
for i, gamma in enumerate(t):
  hself[i] = Hex.hard_moduli_asaro1983(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop)[nactive]
  hlatnt[:,i] = Hex.hard_moduli_asaro1983(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop)
  dhself[i] = Hex.hard_moduli_asaro1983_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop)[nactive]
  dhlatnt[:,i] = Hex.hard_moduli_asaro1983_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop)

  hself1[i] = Hex.hard_moduli_Bassani1991(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop1)[nactive]
  hlatn[:,i] = Hex.hard_moduli_Bassani1991(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop1)
  dhself1[:,i] = Hex.hard_moduli_Bassani1991_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop1)[nactive]
  dhlatn[:,:,i] = Hex.hard_moduli_Bassani1991_dhself(gamma, gmsltl, gamtol, slipnor, slipdir, nactive, prop1)

  gmsltl[nactive] = gamma
  gamtol = np.sum(np.abs(gmsltl))
print(slipnor[:,nactive], slipdir[:,nactive])
# plotting slip system
fig, ax = plt.subplots(1, 1, subplot_kw={"projection": "3d"})
ax.set_proj_type("ortho")
ax.scatter(csys[0], csys[1], csys[2], marker='.',c='b')
ax.scatter(slipdir[0], slipdir[1], zs=slipdir[2], marker='x',c='r')
ax.scatter(slipnor[0], slipnor[1], zs=slipnor[2], marker='*',c='g')
plt.show()

# plotting hardening function
fig, ax = plt.subplots(1, 1)
#ax.plot(t, hself, label='self hardening')
#ax.plot(t, np.transpose(np.delete(hlatnt, nactive, axis=0)), label='latent hardening')
#ax.plot(t, dhself, label='dself')
#ax.plot(t, np.transpose(np.delete(dhlatnt, nactive, axis=0)), label='dlatent')

#ax.plot(t, hself1, label='self hardening')
#ax.plot(t, dhself1[nactive], label='dself')
#ax.plot(t, np.transpose(np.delete(hlatn, nactive, axis=0)), label='latent hardening')
#ax.plot(t, np.transpose(np.delete(dhlatn[:,nactive,:], nactive, axis=0)), label='dlatent')
for j in range(len(gmsltl)):
  ax.plot(t, np.transpose(dhself1[j]), label=j+1, )
#ax.plot(t, np.transpose(dhlatn[nactive,:,:]), label='dlatent')
ax.legend()
#ax.grid(True)
plt.show()

#fig = plt.figure(num=len(t), figsize=(6, 7), dpi=100, facecolor=None, edgecolor=None, frameon=True, clear=False)
#ax.plot(t, np.transpose(dhlatn[nactive,:,:]), label='dlatent')
#im = plt.imshow(dhlatn[:,:,0], cmap='jet')
#print(dhlatn[:, nactive, 0])
#kwargs = {'norm': colors.SymLogNorm(linthresh=0.03, linscale=10), 'cmap': 'Spectral', 'interpolation': None,}
kwargs = {'cmap': 'jet', 'interpolation': None,}
im = plt.matshow(dhlatn[:, :, 0], fignum=0, **kwargs)
cbar = plt.colorbar(im, shrink=0.7)
cbar.formatter.set_powerlimits((0, 0))
im.axes.set_xlabel('$dh_{self}$')
im.axes.set_ylabel('$dh_{latn}$')
plt.show()

#writer = ImageMagickWriter(fps=20)
#fig = plt.figure()
#with writer.saving(fig, "dh.tiff", 100):
#  for j in range(len(t)):
#    im = plt.matshow(dhlatn[:, :, j], fignum=j, **kwargs)
#    writer.grab_frame()

