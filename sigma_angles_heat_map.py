import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import seaborn as sns

num_frames = 1 # number of frames on which perform the analysis
n_chains = 4
n_beads_chain = 10
n_beads_tot = n_chains*n_beads_chain # total number of beads in the system
frame_hm = 0

def distance(x1,y1,z1,x2,y2,z2): # distance between two subsequent beads
    return np.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))

def cosine(dx1,dy1,dz1,dx2,dy2,dz2): # bending potential
    cos = np.dot([dx1,dy1,dz1],[dx2,dy2,dz2])/(np.linalg.norm([dx1,dy1,dz1])*np.linalg.norm([dx2,dy2,dz2]))
    return cos

initial_frame_line = list(range(9, 2000*num_frames, n_beads_tot+9)) # initial lines in the dump_no_bar file from which retrieve the coordinates of the beads
initial_frame_lines_box = list(range(6, 75*num_frames, n_beads_tot+9)) # initial lines in the dump_no_bar file from which retrieve the coordinates of the box

time_steps = np.zeros(num_frames)
cos_frame = np.zeros((n_chains,n_beads_chain))
R0 = 1.5

# the matrix is defined for the whole trajectory
box_cols = 2
box_rows = 3*num_frames
box_matrix = np.zeros((box_rows, box_cols)) 

# for the whole trajectory the length_matrix is defined, with the length of the sides of the box 
length_cols = 3 # Lx,Ly,Lz
length_rows = num_frames
length_matrix = np.zeros((length_rows, length_cols))

j=0

with open('dump_no_bar.lammpstrj', 'r') as file: # to initialize the box_matrix with the coordinates of the sides of the box
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
with open('dump_no_bar.lammpstrj', 'r') as file: # calculate the distance matrix and the bonded energies
    lines = file.readlines()

    for i in initial_frame_line[0:num_frames]: 
        time_steps[k]=np.array(lines[i-8],dtype=float)
        frame = lines[i:i+n_beads_tot]
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float) # In data_T there are all the coordinates and the identification numbers for each atom for the whole trajectory
        sorted_indices = np.argsort(data_T[0,:])
        data_T = data_T[:,sorted_indices] # to order the coordinates from 1 to 40
        
        dist = np.zeros((n_beads_tot, 3)) # distance along x,y,z between subsequent beads in all the chains for a single frame
                
        for n in range(n_chains):
            if n == 0: 
                for j in range(9):
                    if j==2:
                        dist[j][0] = data_T[2][j+1]-data_T[2][j]
                        dist[j][1] = data_T[3][j+1]-data_T[3][j]+length_matrix[k][1]
                        dist[j][2] = data_T[4][j+1]-data_T[4][j]
                    elif j==5:
                        dist[j][0] = data_T[2][j+1]-data_T[2][j]+length_matrix[k][0]
                        dist[j][1] = data_T[3][j+1]-data_T[3][j]
                        dist[j][2] = data_T[4][j+1]-data_T[4][j]
                    else: 
                        dist[j][0] = data_T[2][j+1]-data_T[2][j]
                        dist[j][1] = data_T[3][j+1]-data_T[3][j]
                        dist[j][2] = data_T[4][j+1]-data_T[4][j]

                        if distance(data_T[2][j],data_T[3][j],data_T[4][j],data_T[2][j+1],data_T[3][j+1],data_T[4][j+1])>R0:
                            dist[j][2] = data_T[4][j+1]-data_T[4][j]-length_matrix[k][2]
                        
                dist[9][0] = data_T[2][9]-data_T[2][0]
                dist[9][1] = data_T[3][9]-data_T[3][0]
                dist[9][2] = data_T[4][9]-data_T[4][0]
                if distance(data_T[2][9],data_T[3][9],data_T[4][9],data_T[2][0],data_T[3][0],data_T[4][0])>R0:
                    dist[9][2] = data_T[4][9]-data_T[4][0]
            elif n == 1: 
                for j in range(9):
                    if j==1:
                        dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]
                        dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]-length_matrix[k][1]
                        dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]
                    else:
                        dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]
                        dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]
                        dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]

                        if distance(data_T[2][j+10*n],data_T[3][j+10*n],data_T[4][j+10*n],data_T[2][j+1+10*n],data_T[3][j+1+10*n],data_T[4][j+1+10*n])>R0:
                            dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]-length_matrix[k][2]
  
                dist[9+10*n][0] = data_T[2][19]-data_T[2][10]-length_matrix[k][0]
                dist[9+10*n][1] = data_T[3][19]-data_T[3][10]
                dist[9+10*n][2] = data_T[4][19]-data_T[4][10]+length_matrix[k][2]

                if distance(data_T[2][19],data_T[3][19],data_T[4][19],data_T[2][10]-length_matrix[k][0],data_T[3][10],data_T[4][10]+length_matrix[k][2])  >R0:
                    dist[9+10*n][2] = data_T[4][19]-data_T[4][10]

            elif n == 2: 
                for j in range(9):
                    if j ==7:
                       
                        dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]+length_matrix[k][0]
                        dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]
                        dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]
                    else:
                        dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]
                        dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]
                        dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]

                        if distance(data_T[2][j+10*n],data_T[3][j+10*n],data_T[4][j+10*n],data_T[2][j+1+10*n],data_T[3][j+1+10*n],data_T[4][j+1+10*n]) >R0:
                            dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]+length_matrix[k][2]

                dist[9+10*n][0] = data_T[2][29]-data_T[2][20]
                dist[9+10*n][1] = data_T[3][29]-data_T[3][20]-length_matrix[k][1]
                dist[9+10*n][2] = data_T[4][29]-data_T[4][20]+length_matrix[k][2]

                if distance(data_T[2][29],data_T[3][29],data_T[4][29],data_T[2][20],data_T[3][20]-length_matrix[k][1],data_T[4][20]+length_matrix[k][2])>R0:
                    dist[9+10*n][2] = data_T[4][29]-data_T[4][20]

            else: 
                for j in range(9):
                    dist[j+10*n][0] = data_T[2][j+1+10*n]-data_T[2][j+10*n]
                    dist[j+10*n][1] = data_T[3][j+1+10*n]-data_T[3][j+10*n]
                    dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]
                    if distance(data_T[2][j+10*n],data_T[3][j+10*n],data_T[4][j+10*n],data_T[2][j+1+10*n],data_T[3][j+1+10*n],data_T[4][j+1+10*n])>R0:
                        dist[j+10*n][2] = data_T[4][j+1+10*n]-data_T[4][j+10*n]+length_matrix[k][2]

                dist[10*n][0] = data_T[2][31]-data_T[2][30]-length_matrix[k][0]
                dist[10*n][1] = data_T[3][31]-data_T[3][30]
                dist[10*n][2] = data_T[4][31]-data_T[4][30] 

                dist[8+10*n][0] = data_T[2][39]-data_T[2][38]
                dist[8+10*n][1] = data_T[3][39]-data_T[3][38]
                dist[8+10*n][2] = data_T[4][39]-data_T[4][38]+length_matrix[k][2]
                if distance(data_T[2][38],data_T[3][38],data_T[4][38],data_T[2][39],data_T[3][39],data_T[4][39]+length_matrix[k][2]) >R0:
                    dist[8+10*n][2] = data_T[4][39]-data_T[4][38]

                dist[9+10*n][0] = data_T[2][39]-data_T[2][30]
                dist[9+10*n][1] = data_T[3][39]-data_T[3][30]+length_matrix[k][1]
                dist[9+10*n][2] = data_T[4][39]-data_T[4][30]

                for n in range(n_chains):
                    for m in range(n_beads_chain-1):
                        cos_frame[n][m] = cosine(dist[m+10*n][0],dist[m+10*n][1],dist[m+10*n][2],dist[m+1+10*n][0],dist[m+1+10*n][1],dist[m+1+10*n][2])
                    cos_frame[n][9] = cosine(dist[9+10*n][0],dist[9+10*n][1],dist[9+10*n][2],dist[10*n][0],dist[10*n][1],dist[10*n][2])
        print(f"cosine matrix: {cos_frame}")
        if k == frame_hm:
            ax = sns.heatmap(cos_frame, linewidth=0.5)
            plt.show()

        k+=1 # switch to the following frame
