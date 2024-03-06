import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n = 32
with open('Energy_072.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    E_b=np.zeros(n)
    E_nb=np.zeros(n)

    for k in range(n):
        E_b[k] = float(data_T[0][k])
        E_nb[k] = float(data_T[1][k])

with open('Energy_073.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    E_b1=np.zeros(n)
    E_nb1=np.zeros(n)

    for k in range(n):
        E_b1[k] = float(data_T[0][k])
        E_nb1[k] = float(data_T[1][k])

with open('Energy_074.data', 'r') as file:
    lines = file.readlines()
    data = (l.strip().split(',') for l in lines)
    data_T = np.array(list(zip(*data)),dtype=float)
    E_b2=np.zeros(n)
    E_nb2=np.zeros(n)

    for k in range(n):
        E_b2[k] = float(data_T[0][k])
        E_nb2[k] = float(data_T[1][k])


x = np.zeros(n)
for k in range(n):
    x[k] = 1+(k)*0.01

fig, ax = plt.subplots() # plot with the length of the box along y and z
ax.plot(x, E_b, color = 'm', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.72')
ax.plot(x, E_b1, color = 'c', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.73')
ax.plot(x, E_b2, color = 'orange', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.74')
ax.legend()
ax.set_title("Bonded interaction", fontdict = dict(fontsize = 20))
ax.set_xlabel("Elongation", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$E_{b}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.05))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.01))

plt.show()

fig, ax = plt.subplots() # plot with the length of the box along y and z
ax.plot(x, E_nb, color = 'm', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.72')
ax.plot(x, E_nb1, color = 'c', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.73')
ax.plot(x, E_nb2, color = 'orange', marker = '.',markersize=0.001, label = r'$\sigma_{wca}$ = 0.74')
ax.legend()
ax.set_title("Non bonded interaction", fontdict = dict(fontsize = 20))
ax.set_xlabel("Elongation", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$E_{nb}$",fontdict = dict(fontsize = 16))
plt.grid(visible=None, which='major', axis='both')

ax.xaxis.set_major_locator(MultipleLocator(.05))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(.01))

plt.show()