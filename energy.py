import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from numpy.linalg import norm

num_frames = 90 # number of frames on which perform the analysis
n_chains = 6
n_beads_chain = 10
n_beads_tot = n_chains*n_beads_chain # total number of beads in the system

def distance(x1,y1,z1,x2,y2,z2): # distance between two subsequent beads
    return np.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))

def fene(r,K,R0,epsilon,sigma): # fene potential
    return -0.5*K*(R0**2)*(math.log(1-((r/R0)**2)))+4*epsilon*(((sigma/r)**12)-((sigma/r)**6))+epsilon

def wca(r,epsilon,sigma): # wca potential
    return 4*epsilon*(((sigma/r)**(12))-((sigma/r)**6))+epsilon

def bending(dx1,dy1,dz1,dx2,dy2,dz2,k): # bending potential
    cos = np.dot([dx1,dy1,dz1],[dx2,dy2,dz2])/(np.linalg.norm([dx1,dy1,dz1])*np.linalg.norm([dx2,dy2,dz2]))
    return k*(1-cos)

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
epsilon_wca = 1
sigma_wca = 2.2
rc = sigma_wca*(2**(1/6)) # cutoff
print("rc = ",rc)

E_fene = np.zeros((n_chains,num_frames)) # fene matrix
E_bend = np.zeros((n_chains,num_frames)) # bending matrix
E_b = np.zeros(num_frames) # bonded energy (fene+bending)
E_nb = np.zeros(num_frames) # non bonded energy (wca)
time_steps = np.zeros(num_frames)
# distances between subsequent atoms in the same chain
r = np.zeros((n_chains,n_beads_chain)) 

# the matrix is defined for the whole trajectory
box_cols = 2
box_rows = 3*num_frames
box_matrix = np.zeros((box_rows, box_cols)) 

# for the whole trajectory the length_matrix is defined, with the length of the sides of the box 
length_cols = 3 # Lx,Ly,Lz
length_rows = num_frames
length_matrix = np.zeros((length_rows, length_cols))

# distance between subsequent beads in all the chains for a single frame
dist_rows = n_beads_tot 
dist_cols = 3
dist = np.zeros((dist_rows, dist_cols)) 

j=0

with open('dump.lammpstrj', 'r') as file: # to initialize the box_matrix with the coordinates of the sides of the box
    lines = file.readlines()
    
    for i in initial_frame_lines_box[0:num_frames]: 
        frame = lines[i-1:i+2]
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float)
        sorted_indices = np.argsort(data_T[0,:])
        data_T = data_T[:,sorted_indices] # to order the beads from 1 to 60

        for k in range(2):
            box_matrix[j][k] = float(data_T[k][0])
            box_matrix[j+1][k] = float(data_T[k][1])
            box_matrix[j+2][k] = float(data_T[k][2])

        j=j+3

#print("box matrix: \n",box_matrix)

j=0
for k in range(length_rows):
    length_matrix[k][0] = box_matrix[j][1]-box_matrix[j][0]
    length_matrix[k][1] = box_matrix[j+1][1]-box_matrix[j+1][0]
    length_matrix[k][2] = box_matrix[j+2][1]-box_matrix[j+2][0]
    j=j+3

k=0 # to count the frame
with open('dump.lammpstrj', 'r') as file: # calculate the distance matrix and the bonded energies
    lines = file.readlines()

    for i in initial_frame_line[0:num_frames]: 
        time_steps[k]=np.array(lines[i-8],dtype=float)
        frame = lines[i:i+n_beads_tot]
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float) # In data_T there are all the coordinates and the identification numbers for each atom for the whole trajectory
        sorted_indices = np.argsort(data_T[0,:])
        data_T = data_T[:,sorted_indices] # to order the coordinates from 1 to 60
        
        dist = np.zeros((dist_rows, dist_cols)) 
                
        for n in range(n_chains):
            for j in range(8):
                r[n][j+1] = distance(data_T[2][j+10*n],data_T[3][j+10*n],data_T[4][j+10*n],data_T[2][j+1+10*n],data_T[3][j+1+10*n],data_T[4][j+1+10*n])
                dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]
                dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]
                dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]

        for j in range(n_chains): # bonds 10-1 for all chains
            r[j][0] = distance(data_T[2][j*10],data_T[3][j*10],data_T[4][j*10],data_T[2][j*10+9],data_T[3][j*10+9],data_T[4][j*10+9])
            dist[9+10*j][0] = data_T[2][10*j]-data_T[2][9+10*j]
            dist[9+10*j][1] = data_T[3][10*j]-data_T[3][9+10*j]
            dist[9+10*j][2] = data_T[4][10*j]-data_T[4][9+10*j]
        
        for n in range(n_chains):
            if n == 0 or n==1:
                r[n][9] = distance(data_T[2][8+10*n],data_T[3][8+10*n],data_T[4][8+10*n],data_T[2][9+10*n],data_T[3][9+10*n]+length_matrix[k][1],data_T[4][9+10*n]) # chains along y
                dist[8+10*n][0] = data_T[2][8+1+10*n]-data_T[2][8+10*n]
                dist[8+10*n][1] = data_T[3][8+1+10*n]-data_T[3][8+10*n]+length_matrix[k][1]
                dist[8+10*n][2] = data_T[4][8+1+10*n]-data_T[4][8+10*n]
            if n == 2 or n==3:
                r[n][9] = distance(data_T[2][8+10*n],data_T[3][8+10*n],data_T[4][8+10*n],data_T[2][9+10*n],data_T[3][9+10*n],data_T[4][9+10*n]+length_matrix[k][2]) # chains along z
                dist[8+10*n][0] = data_T[2][8+1+10*n]-data_T[2][8+10*n]
                dist[8+10*n][1] = data_T[3][8+1+10*n]-data_T[3][8+10*n]
                dist[8+10*n][2] = data_T[4][8+1+10*n]-data_T[4][8+10*n]+length_matrix[k][2]
            if n == 4 or n==5:
                r[n][9] = distance(data_T[2][8+10*n],data_T[3][8+10*n],data_T[4][8+10*n],data_T[2][9+10*n]+length_matrix[k][0],data_T[3][9+10*n],data_T[4][9+10*n]) # chains along x
                dist[8+10*n][0] = data_T[2][8+1+10*n]-data_T[2][8+10*n]+length_matrix[k][0]
                dist[8+10*n][1] = data_T[3][8+1+10*n]-data_T[3][8+10*n]
                dist[8+10*n][2] = data_T[4][8+1+10*n]-data_T[4][8+10*n]       

        #print("distance: \n",dist)

        for m in range(n_chains):
            for n in range(n_beads_chain): # fene energy for the frame k
                E_fene[m][k] = E_fene[m][k] + fene(r[m][n],K,R0,epsilon_fene,sigma_fene)
            for n in range(8):
                E_bend[m][k] = E_bend[m][k] + bending(dist[n+10*m][0],dist[n+10*m][1],dist[n+10*m][2],dist[n+1+10*m][0],dist[n+1+10*m][1],dist[n+1+10*m][2],k_bend)
            E_bend[m][k] = E_bend[m][k] + bending(dist[9+10*m][0],dist[9+10*m][1],dist[9+10*m][2],dist[10*m][0],dist[10*m][1],dist[10*m][2],k_bend)
            
        for j in range(n_chains):
            E_b[k] = E_b[k] + E_fene[j][k] + E_bend[j][k]

        E_wca_frame = np.zeros((n_chains-1,n_chains-1)) # matrix with the non bonded interaction for each couple of chains in the system

        for i in range(n_chains-1):
            for j in range(i+1,n_chains):
                dist_2chain=np.zeros((n_beads_chain,n_beads_chain))
                for m in range(n_beads_chain):
                    for n in range(n_beads_chain):
                        dist_2chain[m][n]=distance(data_T[2][m+10*i],data_T[3][m+10*i],data_T[4][m+10*i],data_T[2][n+10*j],data_T[3][n+10*j],data_T[4][n+10*j])
                        if dist_2chain[m][n]<rc:
                            E_wca_frame[i][j-(i+1)] = E_wca_frame[i][j-(i+1)] + wca(dist_2chain[m][n],epsilon_wca,sigma_wca)
                print(f"distance between beads of chains {i} and {j}: \n",dist_2chain)
        
        E_nb[k] += np.sum(E_wca_frame) # sum the contributions from all the couples of chains in the system
        
        k+=1 # switch to the following frame

print("E_non_bond: \n", E_nb)
#print("time steps: \n", time_steps)

min = 0 # initial frame
max = 92 # final frame

fig, ax = plt.subplots()
ax.scatter(time_steps[min:max], E_b[min:max], color = 'b', marker = '2', label = r'$\Delta$t = 0.001')
ax.legend()
ax.set_title("Bonded Energy", fontdict = dict(fontsize = 20))
ax.set_xlabel(r"timestep [$\tau$]", fontdict = dict(fontsize = 16))
ax.set_ylabel(r"$U_{bonded}$[$\epsilon$]",fontdict = dict(fontsize = 16))
          
ax.xaxis.set_major_locator(MultipleLocator(2000))
ax.xaxis.set_major_formatter('{x:.08}')
ax.xaxis.set_minor_locator(MultipleLocator(200))
#plt.show()
