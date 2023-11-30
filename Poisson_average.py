import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n = 9
with open('L_average.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
        #print(data_T)
    Lx=np.zeros(n)
    Ly=np.zeros(n)
    Lz=np.zeros(n)

    for k in range(n):
        Lx[k] = float(data_T[0][k])
        Ly[k] = float(data_T[1][k])
        Lz[k] = float(data_T[2][k])

x = np.zeros(n)
for k in range(n):
    x[k] = 1+k*0.01

Poisson = np.zeros((n-1,2)) # matrix with the poisson's ratio

for k in range(n-1):
    Poisson[k][0] = -((Ly[k+1]-Ly[k])*Lx[k])/((Lx[k+1]-Lx[k])*Ly[k])
    Poisson[k][1] = -((Lz[k+1]-Lz[k])*Lx[k])/((Lx[k+1]-Lx[k])*Lz[k])

#print(Poisson)

fig, ax = plt.subplots() # plot with the length of the box along y and z
ax.plot(x, Ly, color = 'b', marker = '2', label = r'$L_{y}$')
ax.plot(x, Lz, color = 'm', marker = '2', label = r'$L_{z}$')
ax.legend()
#ax.set_title("", fontdict = dict(fontsize = 20))
ax.set_xlabel("Size of the box along y and z", fontdict = dict(fontsize = 16))
ax.set_ylabel("elongation",fontdict = dict(fontsize = 16))
plt.show()

fig, ax = plt.subplots() # plot with the length of the box along y and z
ax.plot(x[1:n], Poisson[:,0], color = 'b', marker = '2', label = r'$L_{y}$')
ax.plot(x[1:n], Poisson[:,1], color = 'm', marker = '2', label = r'$L_{z}$')
ax.legend()
ax.set_xlabel("Poisson's ratio", fontdict = dict(fontsize = 16))
ax.set_ylabel("elongation",fontdict = dict(fontsize = 16))
plt.show()
