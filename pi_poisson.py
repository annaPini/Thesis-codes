import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n = 11

with open('L_average_24.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    Lx=np.zeros(n)
    Ly=np.zeros(n)
    Lz=np.zeros(n)

    for k in range(n):
        Lx[k] = float(data_T[0][k])
        Ly[k] = float(data_T[1][k])
        Lz[k] = float(data_T[2][k])

Poisson = np.zeros((n-1,2))

for k in range(n-1):
    Poisson[k][0] = -((Ly[k+1]-Ly[k])*Lx[k])/((Lx[k+1]-Lx[k])*Ly[k])
    Poisson[k][1] = -((Lz[k+1]-Lz[k])*Lx[k])/((Lx[k+1]-Lx[k])*Lz[k])

x = np.zeros(n)
for k in range(n):
    x[k] = 1+(k)*0.01

x5 = np.zeros(n)
for k in range(n):
    x5[k] = 1+(k)*0.005

fig, ax = plt.subplots() 

ax.plot(x, Ly, color = 'm', marker = '.',markersize=0.001, label = r'$L_{y}$ for $\sigma_{WCA} = 2.4$$\sigma$')
ax.plot(x, Lz, color = 'blue', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 2.4$$\sigma$')
ax.legend()
ax.set_title(r"Size of the box", fontdict = dict(fontsize = 20))
ax.set_xlabel("Multiplicative factor", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$L_{z}$,$L_{y}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.02))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.005))

plt.show()

fig, ax = plt.subplots()
ax.plot(x[1:n], Poisson[:,0], color = 'm', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 2.4$$\sigma$')
ax.plot(x[1:n], Poisson[:,1], color = 'blue', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 2.4$$\sigma$')
ax.legend()
ax.set_title(r"Poisson's ratio", fontdict = dict(fontsize = 20))
ax.set_xlabel("Multiplicative factor", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$\nu_{zx}$,$\nu_{yx}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.02))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.005))

plt.show()