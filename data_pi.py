import numpy as np

initial_line = 10
n_beads = 60
n_chains = 6
Nx = 1
Ny = 1
Nz = 1
coor_box = np.zeros((3,2))
size_box = np.zeros(3)
bonds = np.zeros((n_beads*Nx*Ny*Nz,2))

with open('dump_no_bar_test.lammpstrj', 'r') as file: 
    lines = file.readlines()

    box = lines[5:8]
    data = (l.strip().split(' ') for l in box)
    data_T = np.array(list(zip(*data)),dtype=float) 
    for i in range(3):
        coor_box[i][0] = data_T[0][i]
        coor_box[i][1] = data_T[1][i]

for i in range(3):
    size_box[i] = coor_box[i][1]-coor_box[i][0]

with open('dump_no_bar_test.lammpstrj', 'r') as file: 
    lines = file.readlines()

    frame = lines[initial_line-1:initial_line+n_beads-1]
    data = (l.strip().split(' ') for l in frame)
    data_T = np.array(list(zip(*data)),dtype=float) 
    sorted_indices = np.argsort(data_T[0,:])
    data_T = data_T[:,sorted_indices]

with open('pi_large.data', 'a') as file:
    file.write(f"# LAMMPS data file\n\n{n_beads*Nx*Ny*Nz} atoms\n{(n_beads-n_chains)*Nx*Ny*Nz} bonds\n{(n_beads-2*n_chains)*Nx*Ny*Nz} angles\n\n1 atom types\n1 bond types\n1 angle types\n\n{coor_box[0][0]} {coor_box[0][1]*Nx*Ny*Nz} xlo xhi\n{coor_box[0][0]} {coor_box[0][1]*Nx*Ny*Nz} xlo xhi\n{coor_box[0][0]} {coor_box[0][1]*Nx*Ny*Nz} xlo xhi\n\nMasses\n\n1 1\n\nAtoms\n\n")
    n=0
    for i in range(Nx):
        for j in range(Ny):
            for l in range(Nz):
                for k in range(n_beads):
                    line = ' '.join(map(str, [int(data_T[0][k])+n_beads*n, int(data_T[1][k]), int(data_T[1][k]),data_T[2][k]+i*size_box[0],data_T[3][k]+j*size_box[1],data_T[4][k]+l*size_box[2]]))
                    file.write(line + "\n")
                n = n+1
    file.write("\n\nBonds\n\n")
    for n in range(Nx*Ny*Nz):
        for j in range(n_chains):
            for k in range(n_beads//n_chains-1):
                line = ' '.join(map(str, [k+1+9*j+(n_beads-n_chains)*n,1,k+1 +10*j+n_beads*n,k+2+10*j+n_beads*n]))
                file.write(line + "\n")
    file.write("\nAngles\n\n")
    for n in range(Nx*Ny*Nz):
        for j in range(n_chains):
            for k in range(n_beads//n_chains-2):
                line = ' '.join(map(str, [k+1+8*j+(n_beads-2*n_chains)*n,1,k+1+10*j+n_beads*n,k+2+10*j+n_beads*n,k+3+10*j+n_beads*n]))
                file.write(line + "\n")
file.close()
