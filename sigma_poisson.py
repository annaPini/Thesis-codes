import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n = 32

with open('L_average_72.data', 'r') as file:
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

with open('L_average_73.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
        #print(data_T)
    Lx1=np.zeros(n)
    Ly1=np.zeros(n)
    Lz1=np.zeros(n)

    for k in range(n):
        Lx1[k] = float(data_T[0][k])
        Ly1[k] = float(data_T[1][k])
        Lz1[k] = float(data_T[2][k])

Poisson1 = np.zeros((n-1,2))
for k in range(n-1):
    Poisson1[k][0] = -((Ly1[k+1]-Ly1[k])*Lx1[k])/((Lx1[k+1]-Lx1[k])*Ly1[k])
    Poisson1[k][1] = -((Lz1[k+1]-Lz1[k])*Lx1[k])/((Lx1[k+1]-Lx1[k])*Lz1[k])

with open('L_average_74.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    Lx2=np.zeros(n)
    Ly2=np.zeros(n)
    Lz2=np.zeros(n)

    for k in range(n):
        Lx2[k] = float(data_T[0][k])
        Ly2[k] = float(data_T[1][k])
        Lz2[k] = float(data_T[2][k])

Poisson2 = np.zeros((n-1,2)) 

for k in range(n-1):
    Poisson2[k][0] = -((Ly2[k+1]-Ly2[k])*Lx2[k])/((Lx2[k+1]-Lx2[k])*Ly2[k])
    Poisson2[k][1] = -((Lz2[k+1]-Lz2[k])*Lx2[k])/((Lx2[k+1]-Lx2[k])*Lz2[k])

with open('L_average_75.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    Lx3=np.zeros(n)
    Ly3=np.zeros(n)
    Lz3=np.zeros(n)

    for k in range(n):
        Lx3[k] = float(data_T[0][k])
        Ly3[k] = float(data_T[1][k])
        Lz3[k] = float(data_T[2][k])

Poisson3 = np.zeros((n-1,2)) # matrix with the poisson's ratio

for k in range(n-1):
    Poisson3[k][0] = -((Ly3[k+1]-Ly3[k])*Lx3[k])/((Lx3[k+1]-Lx3[k])*Ly3[k])
    Poisson3[k][1] = -((Lz3[k+1]-Lz3[k])*Lx3[k])/((Lx3[k+1]-Lx3[k])*Lz3[k])

with open('L_average_85.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    Lx4=np.zeros(n)
    Ly4=np.zeros(n)
    Lz4=np.zeros(n)

    for k in range(n):
        Lx4[k] = float(data_T[0][k])
        Ly4[k] = float(data_T[1][k])
        Lz4[k] = float(data_T[2][k])

Poisson4 = np.zeros((n-1,2))

for k in range(n-1):
    Poisson4[k][0] = -((Ly4[k+1]-Ly4[k])*Lx4[k])/((Lx4[k+1]-Lx4[k])*Ly4[k])
    Poisson4[k][1] = -((Lz4[k+1]-Lz4[k])*Lx4[k])/((Lx4[k+1]-Lx4[k])*Lz4[k])

x = np.zeros(n)
for k in range(n):
    x[k] = 1+(k)*0.01

x5 = np.zeros(n)
for k in range(n):
    x5[k] = 1+(k)*0.005

fig, ax = plt.subplots() 

#ax.plot(x5, Ly, color = 'b', marker = '.',markersize=0.001, label = r'$L_{y}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x, Lz, color = 'blue', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 0.72$')
#ax.plot(x5, Ly, color = 'b', marker = '.',markersize=0.001, label = r'$L_{y}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x, Lz1, color = 'r', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 0.73$')
#ax.plot(x5, Ly, color = 'b', marker = '.',markersize=0.001, label = r'$L_{y}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x, Lz2, color = 'm', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 0.74$')
#ax.plot(x, Ly1, color = 'g', marker = '.',markersize=0.001, label = r'$L_{y}$ for $\sigma_{WCA} = 0.8$')
ax.plot(x, Lz3, color = 'c', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 0.75$')
#ax.plot(x, Ly2, color = 'pink', marker = '.', markersize=0.001,label = r'$L_{y}$ for $\sigma_{WCA} = 0.85$')
ax.plot(x, Lz4, color = 'orange', marker = '.',markersize=0.001, label = r'$L_{z}$ for $\sigma_{WCA} = 0.85$')
ax.legend()
ax.set_title(r"Size of the box $L_{z}$", fontdict = dict(fontsize = 20))
ax.set_xlabel("Multiplicative factor", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$L_{z}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.02))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.005))

plt.show()

fig, ax = plt.subplots()
#ax.plot(x5[1:n], Poisson[:,0], color = 'b', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x[1:n], Poisson[:,1], color = 'blue', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 0.72$')
#ax.plot(x5[1:n], Poisson[:,0], color = 'b', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x[1:n], Poisson1[:,1], color = 'r', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 0.73$')
#ax.plot(x5[1:n], Poisson[:,0], color = 'b', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 0.75$')
ax.plot(x[1:n], Poisson2[:,1], color = 'm', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 0.74$')
#ax.plot(x[1:n], Poisson1[:,0], color = 'g', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 0.8$')
ax.plot(x[1:n], Poisson3[:,1], color = 'c', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 0.75$')
#ax.plot(x[1:n], Poisson2[:,0], color = 'purple', marker = '.',markersize=0.001, label = r'$\nu_{yx}$ for $\sigma_{WCA} = 0.85$')
ax.plot(x[1:n], Poisson4[:,1], color = 'orange', marker = '.',markersize=0.001, label = r'$\nu_{zx}$ for $\sigma_{WCA} = 0.85$')
ax.legend()
ax.set_title(r"Poisson's ratio", fontdict = dict(fontsize = 20))
ax.set_xlabel("Multiplicative factor", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$\nu_{zx}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.02))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.005))

plt.show()