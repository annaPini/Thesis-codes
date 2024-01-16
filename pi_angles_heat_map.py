in has to be analysed # 2->y 3->z 4->x

def quicksort(arr): 
    if len(arr) <= 1:
        return arrimport matplotlib.pyplot as plt
import numpy as np
import math
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

n_beads_tot = 60 # total number of beads in the system
beads_chain = 10 # number of beads per chain
f1= 11
f2= 25 # Frame to be analysed
num_frames = 1 # number of frames on which perform the analysis
an_chain = 4 # which cha
    else:
        pivot = arr[0]
        less = [x for x in arr[1:] if x <= pivot]
        greater = [x for x in arr[1:] if x > pivot]
        return quicksort(less) + [pivot] + quicksort(greater)


initial_frame_line = list(range(9, 2000*num_frames, n_beads_tot+9))
print("vector of initial lines:\n",initial_frame_line)

iline_f1 = initial_frame_line[f1]
iline_f2 = initial_frame_line[f2]

cos_cols = num_frames
cos_rows = beads_chain
cos_matrix1 = np.zeros((cos_rows, cos_cols))
cos_matrix2 = np.zeros((cos_rows, cos_cols))

k=0
t=0

with open('dump_elong.lammpstrj', 'r') as file:
    lines = file.readlines()

    frame = lines[iline_f1:iline_f1+n_beads_tot]
    data = (l.strip().split(' ') for l in frame)
    data_T = np.array(list(zip(*data)),dtype=float)
    
    x=[]
    y=[]
    z=[]
    for j in range(n_beads_tot):
        if (data_T[1][j]==an_chain):
            x_coor = float(data_T[2][j])
            x.append(x_coor)
            y_coor = float(data_T[3][j])
            y.append(y_coor)
            z_coor = float(data_T[4][j])
            z.append(z_coor)

    quicksort(x)
    quicksort(y)
    quicksort(z)
    # for 2
    #x[0],x[1]=x[1],x[0]
    #y[0],y[1]=y[1],y[0]
    #z[0],z[1]=z[1],z[0]
    

    dist_rows = 3 # x,y,z
    dist_cols = beads_chain-1 # 9 segments for a single chain

    dist = np.zeros((dist_rows, dist_cols))
    for j in range(dist_cols):
        dist[0][j] = x[j+1]-x[j]
    for j in range(dist_cols):
        dist[1][j] = y[j+1]-y[j]
    for j in range(dist_cols):
        dist[2][j] = z[j+1]-z[j]
    
    for j in range(dist_cols):
        norm_a = np.linalg.norm([row[j-1] for row in dist])
        norm_b = np.linalg.norm([row[j] for row in dist])
        cdot= np.dot([row[j-1] for row in dist],[row[j] for row in dist])/(norm_a*norm_b)
        cos_matrix1[j][k] = cdot

    norm_a = np.linalg.norm([row[0] for row in dist])
    norm_b = np.linalg.norm([row[8] for row in dist])
    cdot_if= np.dot([row[0] for row in dist],[row[8] for row in dist])/(norm_a*norm_b)
    cos_matrix1[dist_cols][k] = cdot_if # angle which goes through the boundary condition

    k+=1
        
print("cosine matrix 1 is: \n",cos_matrix1)  

with open('dump_elong.lammpstrj', 'r') as file:
    lines = file.readlines()

    frame = lines[iline_f2:iline_f2+n_beads_tot]
    data = (l.strip().split(' ') for l in frame)
    data_T = np.array(list(zip(*data)),dtype=float)
    
    x=[]
    y=[]
    z=[]
    for j in range(n_beads_tot):
        if (data_T[1][j]==an_chain): 
            x_coor = float(data_T[2][j])
            x.append(x_coor)
            y_coor = float(data_T[3][j])
            y.append(y_coor)
            z_coor = float(data_T[4][j])
            z.append(z_coor)

    quicksort(x)
    quicksort(y)
    quicksort(z)

    dist_rows = 3 # x,y,z
    dist_cols = beads_chain-1 # 9 segments for a single chain

    dist = np.zeros((dist_rows, dist_cols))
    for j in range(dist_cols):
        dist[0][j] = x[j+1]-x[j]
    for j in range(dist_cols):
        dist[1][j] = y[j+1]-y[j]
    for j in range(dist_cols):
        dist[2][j] = z[j+1]-z[j]
    
    for j in range(dist_cols):
        norm_a = np.linalg.norm([row[j-1] for row in dist])
        norm_b = np.linalg.norm([row[j] for row in dist])
        cdot= np.dot([row[j-1] for row in dist],[row[j] for row in dist])/(norm_a*norm_b)
        cos_matrix2[j][t] = cdot

    norm_a = np.linalg.norm([row[0] for row in dist])
    norm_b = np.linalg.norm([row[8] for row in dist])
    cdot_if= np.dot([row[0] for row in dist],[row[8] for row in dist])/(norm_a*norm_b)
    cos_matrix2[dist_cols][t] = cdot_if # angle which goes through the boundary condition

    t+=1
        
print("cosine matrix 2 is: \n",cos_matrix2)  

def cos_to_degree(cos): # converting the cosines in angles
    if cos>0:
        return math.degrees(math.acos(cos))
    else:
        return 180-math.degrees(math.acos(cos))
    
angle_matrix1 = np.array([[cos_to_degree(element) for element in row] for row in cos_matrix1])
angle_matrix2 = np.array([[cos_to_degree(element) for element in row] for row in cos_matrix2])
print(f"Angle matrix:\n{angle_matrix1}")
print(f"Angle matrix:\n{angle_matrix2}")

n_bin_ang = 30 # how many bins in the hystogram with the angles
bin_ang1 = np.zeros((n_bin_ang)) 
bin_ang2 = np.zeros((n_bin_ang)) 

for i in range(cos_rows): # assign each element of the angine matrix to a bin
    for j in range(cos_cols):
        for k in range(n_bin_ang):
            if (k*90/ n_bin_ang) < angle_matrix1[i][j] <= ((k + 1)*90/ n_bin_ang):
                bin_ang1[k]+=1

for i in range(cos_rows): # assign each element of the angine matrix to a bin
    for j in range(cos_cols):
        for k in range(n_bin_ang):
            if (k*90/ n_bin_ang) < angle_matrix2[i][j] <= ((k + 1)*90/ n_bin_ang):
                bin_ang2[k]+=1

norm_bin_ang1 = (1/(cos_cols*cos_rows))*bin_ang1
norm_bin_ang2 = (1/(cos_cols*cos_rows))*bin_ang2

c_bin_ang = np.zeros(n_bin_ang) #central value of each bin (the same for both the frames)
for i in range(n_bin_ang):
    c_bin_ang[i] = (i+1)*90/n_bin_ang
print(f"Centers of the bins in the angle histogram:\n{c_bin_ang}\nsize: {len(c_bin_ang)}")

plt.bar(c_bin_ang,norm_bin_ang1,alpha = 0.5,width=90/n_bin_ang,color='turquoise',edgecolor='mediumblue',linewidth=2,label = f'Frame {f1}') # bars for the first frame
plt.bar(c_bin_ang,norm_bin_ang2,alpha=0.5,width=90/n_bin_ang,color='violet',edgecolor='darkmagenta',linewidth=2, label = f'Frame {f2}') # bars for the second frame
plt.xticks(c_bin_ang)
plt.xlim(0,92)
plt.legend()
plt.title("Angle distribution for a single chain")
plt.xlabel(r"Angles [$\theta$]")
plt.ylabel("Occurrences")
plt.show()
