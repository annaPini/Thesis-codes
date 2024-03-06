import matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n_beads_tot = 40 
num_frames = 109

initial_frame_line = list(range(6, 750*num_frames, n_beads_tot+9))

box_cols = 2
box_rows = 3*num_frames
box_matrix = np.zeros((box_rows, box_cols)) # x_min/max, y_min/max, z_min/max per ogni frame

length_cols = 3
length_rows = num_frames
length_matrix = np.zeros((length_rows, length_cols)) # length of the sides of the box

j=0

with open('dump_no_bar.lammpstrj', 'r') as file:
    lines = file.readlines()
    
    for i in initial_frame_line[0:num_frames]: 
        frame = lines[i-1:i+2]
        data = (l.strip().split(' ') for l in frame)
        data_T = np.array(list(zip(*data)),dtype=float)
        
        for k in range(2):
            box_matrix[j][k] = float(data_T[k][0])
            box_matrix[j+1][k] = float(data_T[k][1])
            box_matrix[j+2][k] = float(data_T[k][2])

        j=j+3

j=0
for k in range(length_rows):
    length_matrix[k][0] = box_matrix[j][1]-box_matrix[j][0] # length of the side along x for frame k
    length_matrix[k][1] = box_matrix[j+1][1]-box_matrix[j+1][0] # length of the side along y for frame k
    length_matrix[k][2] = box_matrix[j+2][1]-box_matrix[j+2][0] # length of the side along z for frame k
    j=j+3
        
n=0
for k in range(num_frames-1):
    if(length_matrix[k+1][2]-length_matrix[k][2]!=0):
        n=n+1 # to count the frames in which the box relaxes along y and z
print(f"n={n}")

Ly_average = 0
Lz_average = 0
for k in range(num_frames-n,num_frames):
    Ly_average = Ly_average + length_matrix[k][1]/n
    Lz_average = Lz_average + length_matrix[k][2]/n

print(f"Ly_av = {Ly_average} Lz_av = {Lz_average}")

with open('L_average_87.data', 'a') as file:
    #file.write(str([length_matrix[num_frames-n][0],Ly_average,Lz_average])+"\n")
    file.write(str([length_matrix[num_frames-1][0],length_matrix[num_frames-1][1],length_matrix[num_frames-1][2]])+"\n")
    file.close()