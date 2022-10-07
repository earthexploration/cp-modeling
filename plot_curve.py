# -*- coding: utf-8 -*-
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import math
import os

filename = "Hwang1d.dat"
data = np.loadtxt(filename)
# print(data)
x = [t[0] for t in data]
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
plt.xlabel(r"$x$")
plt.ylabel(r"$y$")
plt.title("plot y(x)")
plt.tight_layout()
plt.show()
#plt.savefig(fname+".pdf", format="pdf")
