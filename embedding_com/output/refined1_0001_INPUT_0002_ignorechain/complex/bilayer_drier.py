#!/usr/bin/env python

import sys, re, numpy
from numpy.linalg import norm

in_file = sys.argv[1] # GROMACS coordinate file
ifile	= open(in_file, 'r')
lines	= ifile.readlines()
ifile.close()

# Extract Lipid Phosphate (P8) atomic coordinates: Z-axis
P8_z = [float(l.split()[-1]) for l in lines if re.search(r"P8", l)]
P8_z_mean  = numpy.mean(P8_z)
# Average Z-axis postions, upper leaflet
P8_z_Upper = [z for z in P8_z if z>P8_z_mean] 
P8_z_Upper_mean  = numpy.mean(P8_z_Upper)
# Average Z-axis postions, lower leaflet
P8_z_Lower = [z for z in P8_z if z<P8_z_mean]
P8_z_Lower_mean  = numpy.mean(P8_z_Lower)
del P8_z

# Extract Protein atomic coordinates, to calculate centre of mass (COM)
Protein_xyz = []
for l in lines[2:]:
	if re.search(r"POPC", l):
		break
	else:
		Protein_xyz.append(map(float,l.split()[-3:]))
# Work out COM of Protein 
Protein_COM = numpy.mean(numpy.asarray(Protein_xyz).T, axis=1)
Protein_COM = numpy.asarray(Protein_COM)
# Work out radius of smallest sphere centre at COM, without overlapping
dr = 0.1
r  = 0
clashes = []
while len(clashes) == 0:
	r += dr
	clashes = [xyz for xyz in Protein_xyz if norm(numpy.array(xyz)-Protein_COM)<=r] 
R_sphere = r - dr # Smallest radius
del Protein_xyz

# Extract lines up to first SOL atom
first_OW = next(obj for obj in enumerate(lines) if re.search(r"OW",obj[1]))
head	= lines[2:first_OW[0]] # Exclude headline and number of atoms
# Extract box dimensions
box_dim = lines[-1]
# Extract SOL lines and OW indexes
SOL = lines[first_OW[0]:-1]
del lines

# Select Water coordinates outside bilayer, but within pore
N = len(SOL)
SOL_filtered = []
for i in range(0,N,3):
	z = float(SOL[i].split()[-1])
	if (z < P8_z_Lower_mean or z > P8_z_Upper_mean):
		SOL_filtered.append(SOL[i])
		SOL_filtered.append(SOL[i+1])
		SOL_filtered.append(SOL[i+2])
	else:
		r_xy = numpy.array(map(float,SOL[i].split()[-3:-1]))
		if norm(r_xy - Protein_COM[:2]) < R_sphere:
			SOL_filtered.append(SOL[i])
                	SOL_filtered.append(SOL[i+1])
                	SOL_filtered.append(SOL[i+2])

# Print out data into new file
out_file = in_file.split('.')[0]+'_dry.gro'
ofile    = open(out_file, 'w')
N_atoms = str(len(head)+len(SOL_filtered))+'\n'
ofile.write("Protein in dry bilayer in Water\n")
ofile.write(N_atoms)
for l in head:
	ofile.write(l)
for l in SOL_filtered:
	ofile.write(l)
ofile.write(box_dim)
ofile.close()

# Summary of process information
N_waters0 = int(len(SOL)/3)
N_waters = int(len(SOL_filtered)/3)
print("Original number of waters: \t%s"% N_waters0)
print("Waters saved after filtering:\t%d"% N_waters)

