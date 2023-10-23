import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n_beads_tot = 60 # total number of beads in the system
num_frames = 12 # number of frames on which perform the analysis

initial_frame_line = list(range(6, 75*num_frames, n_beads_tot+9))
#print("vector of initial lines:\n",initial_frame_line)
#print("initial line: ",initial_frame_line[0:num_frames])

box_cols = 2
box_rows = 3*num_frames
box_matrix = np.zeros((box_rows, box_cols))

j=0

with open('dump.lammpstrj', 'r') as file:
    lines = file.readlines()
    
    for i in initial_frame_line[0:num_frames]: 
        frame = lines[i-1:i+2]
        #print(frame)
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float)
        #print(data_T)
        
        for k in range(2):
            box_matrix[j][k] = float(data_T[k][0])
            box_matrix[j+1][k] = float(data_T[k][1])
            box_matrix[j+2][k] = float(data_T[k][2])

        j=j+3

length_cols = 3
length_rows = num_frames
length_matrix = np.zeros((length_rows, length_cols))

j=0
for k in range(length_rows):
    length_matrix[k][0] = box_matrix[j][1]-box_matrix[j][0]
    length_matrix[k][1] = box_matrix[j+1][1]-box_matrix[j+1][0]
    length_matrix[k][2] = box_matrix[j+2][1]-box_matrix[j+2][0]
    j=j+3
        
#print("length matrix:\n",length_matrix)

poisson_rows = num_frames-1
poisson_matrix = np.zeros((poisson_rows, 2))
extension_matrix = np.zeros((poisson_rows, 3))


for k in range(poisson_rows):
    ex=(length_matrix[k+1][0]-length_matrix[k][0])
    ey=(length_matrix[k+1][1]-length_matrix[k][1])
    ez=(length_matrix[k+1][2]-length_matrix[k][2])

    if(ex==0):
        extension_matrix[k][0]=0
        extension_matrix[k][1]=0
        extension_matrix[k][2]=0
    else:
        extension_matrix[k][0]= ex
        extension_matrix[k][1] = ey
        extension_matrix[k][2] = ez

print("Extension matrix: \n",extension_matrix)

for k in range(poisson_rows):
    ex=(length_matrix[k+1][0]-length_matrix[k][0])/length_matrix[k][0]
    ey=(length_matrix[k+1][1]-length_matrix[k][1])/length_matrix[k][1]
    ez=(length_matrix[k+1][2]-length_matrix[k][2])/length_matrix[k][2]

    if(ex==0):
        poisson_matrix[k][0]=0
        poisson_matrix[k][1]=0
    else:
        poisson_matrix[k][0]= -ey/ex
        poisson_matrix[k][1] = -ez/ex

#print("Poisson's ratio:\n",poisson_matrix)

x = np.zeros(num_frames-1)

for i in range(num_frames-1):
    x[i]=length_matrix[i][0]-length_matrix[0][0]

fig, ax = plt.subplots()
ax.plot(x, length_matrix[1:,1], color = 'b', marker = '2', label = r'$L_{y}$')
ax.plot(x, length_matrix[1:,2], color = 'm', marker = '2', label = r'$L_{z}$')
ax.legend()
ax.set_title("Lattice Extensions", fontdict = dict(fontsize = 20))
ax.set_xlabel("Elongation along x", fontdict = dict(fontsize = 16))
ax.set_ylabel("Lattice extensions along y and z",fontdict = dict(fontsize = 16))
plt.show()

fig, ax = plt.subplots()
ax.plot(x, poisson_matrix[:,0], color = 'b', marker = '2', label = r'$\nu_{yx}$')
ax.plot(x, poisson_matrix[:,1], color = 'm', marker = '2', label = r'$\nu_{yz}$')
ax.legend()
ax.set_title("Poisson's ratio", fontdict = dict(fontsize = 20))
ax.set_xlabel("Elongation along x", fontdict = dict(fontsize = 16))
ax.set_ylabel("Poisson's ratio",fontdict = dict(fontsize = 16))
plt.show()
