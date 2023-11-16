import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from numpy.linalg import norm

n_beads_tot = 60 # total number of beads in the system
num_frames = 180 # number of frames on which perform the analysis
n_chains = 6
n_beads_chain = 10

def distance(x1,y1,z1,x2,y2,z2): # distance between two beads
    return np.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))

def fene(r): # fene energy
    return -0.5*K*(R0**2)*(math.log(1-(r/R0)**2)) + 4*epsilon_fene*(((sigma_fene/r)**12)-((sigma_fene/r)**6))+epsilon_fene

def bending(x1,y1,z1,x2,y2,z2,k): # bending energy
    a = np.array([x1,y1,z1])
    b = np.array([x2,y2,z2])
    cos = np.dot(a,b)/(norm(a)*norm(b))
    return k*(1+cos)

initial_frame_line = list(range(9, 2000*num_frames, n_beads_tot+9)) # initial lines in the dump file from which retrieve the coordinates of the beads
initial_frame_lines_box = list(range(6, 75*num_frames, n_beads_tot+9)) # initial lines in the dump file from which retrieve the coordinates of the box

# FENE parameters
K=30
R0=1.5
epsilon_fene = 1
sigma_fene = 1
# bending stiffness
k_bend = 20
# WCA parameters
epsilon_wca = 2.2
sigma_wca = 2.4694

E_fene = np.zeros(num_frames) # fene 
E_bend = np.zeros(num_frames) # bending
E_b = np.zeros(num_frames) # bonded energy
E_nb = np.zeros(num_frames) # non bonded energy
r = np.zeros((n_chains,n_beads_chain)) # distances between subsequent atoms in the same chain

box_cols = 2
box_rows = 3*num_frames
box_matrix = np.zeros((box_rows, box_cols))

j=0

with open('dump.lammpstrj', 'r') as file: # to calculate the length of the box for each frame
    lines = file.readlines()
    
    for i in initial_frame_lines_box[0:num_frames]: 
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

k=0
with open('dump.lammpstrj', 'r') as file: # to calculate the distance matrix and the bonded energies
    lines = file.readlines()

    for i in initial_frame_line[0:num_frames]: 
        frame = lines[i:i+n_beads_tot]
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float) # In data_T there are all the coordinates and the identification numbers for each atom for the whole trajectory
        sorted_indices = np.argsort(data_T[0,:])
        data_T = data_T[:,sorted_indices] # to order the coordinates from 1 to 60

        for j in range(n_beads_tot):
            if j>=0 and j<8: # bonds inside the box, from 1-2 to 8-9 for chain 1
                r[0][j+1] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            if j>=10 and j<=17: # bonds inside the box, from 1-2 to 8-9 for chain 2
                r[1][j+1-10] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            if j>=20 and j<=27: # bonds inside the box, from 1-2 to 8-9 for chain 3
                r[2][j+1-20] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            if j>=30 and j<=37: # bonds inside the box, from 1-2 to 8-9 for chain 4
                r[3][j+1-30] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            if j>=40 and j<=47: # bonds inside the box, from 1-2 to 8-9 for chain 5
                r[4][j+1-40] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            if j>=50 and j<=57: # bonds inside the box, from 1-2 to 8-9 for chain 6
                r[5][j+1-50] = distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])
            for j in range(6): # bonds 10-1 for all chains
                r[j][0] = distance(data_T[2][j*10],data_T[3][j*10],data_T[4][j*10],data_T[2][j*10+9],data_T[3][j*10+9],data_T[4][j*10+9])

        # bond 9-10: I need to translate the box
        #r[0][9] = distance(data_T[2][8],data_T[3][8],data_T[4][8],data_T[2][9],data_T[3][9],data_T[4][9])

        E_fene1 = 0 # fene for chain 1
        E_fene2 = 0 # fene for chain 2
        E_fene3 = 0 # fene for chain 3
        E_fene4 = 0 # fene for chain 4
        E_fene5 = 0 # fene for chain 5
        E_fene6 = 0 # fene for chain 6
        E_bend1 = 0 # bending for chain 1
        E_bend2 = 0 # bending for chain 2
        E_bend3 = 0 # bending for chain 3
        E_bend4 = 0 # bending for chain 4
        E_bend5 = 0 # bending for chain 5
        E_bend6 = 0 # bending for chain 6

        for n in range(n_beads_chain-1): # calculation of fene energy for the frame k. The cycle should be in range(n_beads_chain) but it is not bossible to divide for 0 and I did not find how to represent the bond through the pbc
            E_fene1 = E_fene1 + fene(r[0][n])
            E_fene2 = E_fene2 + fene(r[1][n])
            E_fene3 = E_fene3 + fene(r[2][n])
            E_fene4 = E_fene4 + fene(r[3][n])
            E_fene5 = E_fene5 + fene(r[4][n])
            E_fene6 = E_fene6 + fene(r[5][n])
        
        for j in range(n_beads_tot):
            if j>=0 and j<8:
                E_bend1 = E_bend1 + bending(data_T[2][j],data_T[3]  [j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T [4][j+1],k_bend)
            if j>=10 and j<17:
                E_bend2 = E_bend2 + bending(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1],k_bend)
            if j>=20 and j<27:
                E_bend3 = E_bend3 + bending(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1],k_bend)
            if j>=30 and j<37:
                E_bend4 = E_bend4 + bending(data_T[2][j],data_T
                [3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1],k_bend)
            if j>=40 and j<47:
                E_bend5 = E_bend5 + bending(data_T[2][j],data_T
                [3][j],data_T[4][j],data_T[2][j+1],data_T[3][j
                +1],data_T[4][j+1],k_bend)
            if j>=50 and j<57:
                E_bend6 = E_bend6 + bending(data_T[2][j],data_T
                [3][j],data_T[4][j],data_T[2][j+1],data_T[3][j
                +1],data_T[4][j+1],k_bend)

        E_fene[k] = E_fene1 + E_fene2 + E_fene3 + E_fene4 + E_fene5 + E_fene6 # total fene energy in frame k
        E_bend[k] = E_bend1 + E_bend2 + E_bend3 + E_bend4 + E_bend5 + E_bend6
        E_b[k] = E_fene[k]+E_bend[k] 
        k+=1
print("Bonded energy contribution: \n",E_b)
