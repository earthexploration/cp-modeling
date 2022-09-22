# -*- coding: utf-8 -*-
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import os
mpl.use('agg')
filename = "test_asaro2.dat"
data = np.loadtxt(filename)
# print(data)
x = [abs(t[0]) for t in data]
y = [t[1] for t in data]

#sizeOfFont = 12
#fontProperties = {
#'family' : 'serif',
#'serif' : ['Computer Modern Serif'],
#'weight' : 'normal',
#'size' : sizeOfFont
#}
#plt.rc('text', usetex=True)
#plt.rc('font', **fontProperties)

fig = plt.figure()
plt.plot(x,
        y,
        )
plt.grid(True)
plt.xlabel(r"$\gamma$")
plt.ylabel(r"$\tau$")
plt.title("plot y(x)")
plt.tight_layout()
#plt.show()
plt.savefig(filename+".ps", format="ps")
