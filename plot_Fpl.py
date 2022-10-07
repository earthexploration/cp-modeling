# -*- coding: utf-8 -*-
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.figure as Fig
import matplotlib.axes as Ax
import math
import os

def Fpl_c(x, k):
    F = [[1.0/np.cosh(x), np.tanh(x)], [-np.tanh(x), 1.0/np.cosh(x)]]
    y = [k*np.cos(x), k*np.sin(x)]
    z = np.dot(F, y)
    return z

def Fpl_e(x, k, q):
    F = [[1.0/np.cosh(x), np.tanh(x)], [-np.tanh(x), 1.0/np.cosh(x)]]
    r = np.sqrt(k**2*(1-q**2)/(1-q**2*np.cos(x)**2))
    y = [r*np.cos(x), r*np.sin(x)]
    z = np.dot(F, y)
    return z

k = 1
q = 0.9
t = np.linspace(2*math.pi, 4*math.pi, num=1001)
r = np.sqrt(k**2*(1-q**2)/(1-q**2*np.cos(t)**2))
x = []
y = []
for i,tt in enumerate(t): 
    z = Fpl_e(tt, r[i], q)
    #z = Fpl_c(tt, k)
    x.append(z[0])
    y.append(z[1])
#print(x, y)

mpl.use('Qt5Agg')
from pathlib import Path
afm_path = Path(mpl.get_data_path(), 'fonts', 'afm', 'ptmr8a.afm')
from matplotlib.afm import AFM
with afm_path.open('rb') as fh:
    afm = AFM(fh)

sizeOfFont = 12
fontProperties = {
'family' : 'serif',
'serif' : ['Computer Modern Serif'],
'weight' : 'normal',
'size' : sizeOfFont
}
plt.rc('text', usetex=True)
plt.rc('font', **fontProperties)

#fig = Fig.Figure(figsize=[6,8], dpi=200.0, facecolor='black', edgecolor='black', linewidth=0.5, frameon=False, subplotpars=None, tight_layout='tight', constrained_layout='constrained')
#ax = Ax.Axes(fig, [0,0,1,1])
fig = plt.figure()
ax = fig.add_subplot(111)
#line1, = ax.plot(k*np.cos(t), k*np.sin(t), label='initial')
#plt.plot(r*np.cos(t), r*np.sin(t), label='initial')
line1, = ax.plot(r*np.cos(t), r*np.sin(t), label='initial')
line2, = ax.plot(x, y, label='final')

#plt.quiver(X=r*np.cos(t), Y=r*np.sin(t), U=x-r*np.cos(t), V=y-r*np.sin(t))
fig.legend(handles=[line1, line2])
#plt.grid(True)
ax.set_aspect("equal","box")
ax.set_xlabel(r"$x$")
ax.set_ylabel(r"$y$")
ax.set_title("deformation mapping")
#plt.tight_layout()
#fig.savefig(ax.get_title(), transparent=0.1, format='pdf',backend='pgf')
fig.canvas.draw()
fig.canvas.flush_events()
plt.show()

