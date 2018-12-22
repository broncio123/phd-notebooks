#!/bin/bash

path=$1

# Resize simulation box, containing Protein and bilayer atmos only
## Create index file of protein and lipid atomic indexes
if [ ! -f ${path}/protein+lipids.ndx ]; then
	printf "1|13\nq\n"| gmx_mpi make_ndx -f ${path}/confout.gro -o  ${path}/protein+lipids.ndx
fi

## Extract box dimensions
complexfile=${path}/confout.gro
read Lx Ly Lz <<<$(tail -1 $complexfile)
	## Redifine box dimension, Z-axis
Lz_new=10.000
	## Centre bilayer and protein in larger box
printf "24\n24\n"| gmx_mpi editconf -f ${path}/confout.gro -n ${path}/protein+lipids.ndx -o ${path}/protein+lipids.gro -c -box $Lx $Ly $Lz_new

# Resolvate new box, SPC water model
if [ ! -f ${path}/protein+lipids+sol_RAW.gro ]; then 
	gmx_mpi solvate -cp ${path}/protein+lipids.gro -cs spc216.gro -o ${path}/protein+lipids+sol_RAW.pdb
fi
##########################################################################################################################
# Water removal: Outisde box, between lipids in bilayer, and from protein-lipids interface 
## Get PDB of protein alone with chain separation 
if [ ! -f ${path}/protein+lipids+sol_RAW.ndx ]; then
    printf "q\n" | gmx_mpi make_ndx -f ${path}/protein+lipids+sol_RAW.pdb -o ${path}/protein+lipids+sol_RAW.ndx
fi

if [ ! -f ${path}/protein_RAW.pdb ]; then
    echo "Protein" | gmx_mpi editconf -f ${path}/protein+lipids+sol_RAW.pdb -n ${path}/protein+lipids+sol_RAW.ndx -o ${path}/protein_RAW.pdb
    sed -i '/TER/d' ${path}/protein_RAW.pdb
    sed -i 's/ENDMDL/END/g' ${path}/protein_RAW.pdb
    pymol -qc ~/mpmodeling/protocols/sem/label_chains4gro.py -- ${path}/protein_RAW.pdb ${path}/protein_RAW.pdb
fi

##########################################################################################################################

if [ ! -f ${path}/protein+lipids+sol_clean.pdb ]; then
var=${path} python3 <<@@
import os
import numpy, re
import isambard_dev

path = os.environ["var"]

# 1. FIND POSITION OF P8 ATOMS IN POPC
in_file = path+"/"+"protein+lipids+sol_RAW.pdb" 
ifile   = open(in_file, 'r')
lines   = ifile.readlines()
ifile.close()

# Extract Lipid Phosphate (P8) atomic coordinates: Z-axis
P8_z = [float(l.split()[-3]) for l in lines if re.search(r"P8", l)]
P8_z_mean  = numpy.mean(P8_z)
# Average Z-axis postions, upper leaflet
P8_z_Upper = [z for z in P8_z if z>P8_z_mean]
P8_z_Upper_mean  = numpy.mean(P8_z_Upper)
# Average Z-axis postions, lower leaflet
P8_z_Lower = [z for z in P8_z if z<P8_z_mean]
P8_z_Lower_mean  = numpy.mean(P8_z_Lower)
del P8_z

# 2. REMOVE WATERS OUTSIDE BOX
SOL = [l for l in lines if re.search(r'SOL', l)]


# 3. REMOVE WATERS WITHIN BILAYER AND FROM PROTEIN-LIPID INTERFACE 
protein = isambard_dev.ampal.convert_pdb_to_ampal(path+"/"+"protein_RAW.pdb")
protein_com = numpy.array(protein.centre_of_mass)

## Find primitives per chain
ampal = protein
prims = numpy.array([x.coordinates for x in ampal.primitives])
ref_axis = isambard_dev.ampal.pseudo_atoms.Primitive.from_coordinates(numpy.mean(prims, axis=0))

# Determine chain radial profile and define a Z-partition
dist = prims[0] - ref_axis.coordinates
R_profile = numpy.linalg.norm(dist, axis=1)

primitive_z = numpy.array(ref_axis.coordinates).T[2]

protein_z_partition = []
for i in range(len(primitive_z)-1):
    lower_lim = primitive_z[i]
    upper_lim = primitive_z[i+1]
    rad_lim = R_profile[i]
    protein_z_partition.append([lower_lim, upper_lim, rad_lim])

R_low = R_profile[0]
R_upper = R_profile[-1]
total_z_partition = [[P8_z_Lower_mean, primitive_z.min(), R_low]] \
    + protein_z_partition  \
    + [[primitive_z.max(), P8_z_Upper_mean, R_upper]]


## Remove waters inside between bilayer leaflets and from protein-lipid interface,
SOL_data = SOL
N = len(SOL_data)
SOL_inside_indices = []
for s in total_z_partition:
        lower_lim = s[0]
        upper_lim = s[1]
        rad_lim = s[2]
        for i in range(0,N,3):
                x,y,z = list(map(float,SOL_data[i].split()[-5:-2]))
                r_xy = numpy.linalg.norm(numpy.array([x,y]) - protein_com[:2])
                if (lower_lim < z <= upper_lim) and (r_xy < rad_lim):
                        SOL_inside_indices.append(i)

## Find waters outside bilayer,
SOL_data = SOL
N = len(SOL)
SOL_outside_indices = []
for i in range(0,N,3):
        x,y,z = list(map(float,SOL_data[i].split()[-5:-2]))
        r_xy = numpy.linalg.norm(numpy.array([x,y]) - protein_com[:2])
        if (z <= P8_z_Lower_mean) or (z > P8_z_Upper_mean):
                SOL_outside_indices.append(i)

## Merge solvent indeces and select water coordinates accordingly,
SOL_filtered_indeces = SOL_inside_indices + SOL_outside_indices
SOL_filtered_indeces.sort()

# SUMMARY OF WATER REMOVAL:
summary_file = open(path+"/"+"water_removal.txt", 'w')
summary_file.write("Remove waters within bilayer and protein-lipid interface\n")
summary_file.write("Initial number of waters: "+str(len(SOL))+"\n")
summary_file.write("Number of waters within box boundaries: "+str(len(SOL_filtered_indeces))+"\n")
summary_file.write("Water surviving: "+str(len(SOL_filtered_indeces)/float(len(SOL))*100)+" %"+"\n\n")
summary_file.close()

# 4. GENERATE SYSTEM PDB WITH FILTERED WATERS
SOL_filtered = []
for index in SOL_filtered_indeces:
    SOL_filtered.append(SOL[index])
    SOL_filtered.append(SOL[index+1])
    SOL_filtered.append(SOL[index+2])

pdb_prot_popc = []
for l in lines:
        if re.search(r"SOL", l):
                break
        else:
            pdb_prot_popc.append(l)

out_file = path+'/'+'protein+lipids+sol_clean.pdb'
ofile    = open(out_file, 'w')
for l in pdb_prot_popc:
    ofile.write(l)
for l in SOL_filtered:
    ofile.write(l)
ofile.close()
@@
fi

# Centre all atoms in box
gmx_mpi editconf -f ${path}/protein+lipids+sol_clean.pdb -o ${path}/protein+lipids+sol_clean.pdb -c

# Fix topology 
cp ${path}/for_embedding.top ${path}/topol.top
N_filtered=$(grep "Number of waters within box boundaries" ${path}/water_removal.txt | awk -F":" 'END{print $2}')
sed -i "/SOL/d" ${path}/topol.top
sed -i "/K/d" ${path}/topol.top
sed -i "/CL/d" ${path}/topol.top
printf "SOL\t\t $N_filtered\n" >> ${path}/topol.top

# 10. Add K-Cl ions to the system, at 1M concentration, neutral
if [ ! -f ${path}/ionise.tpr ]; then
    gmx_mpi grompp -f ${path}/mdpf/ionise.mdp -c ${path}/protein+lipids+sol_clean.pdb -o ${path}/ionise.tpr -p ${path}/topol.top
fi

if [ ! -f ${path}/ionise.pdb ]; then
    echo "SOL" | gmx_mpi genion -s ${path}/ionise.tpr -p ${path}/topol.top -o ${path}/ionise.gro -conc 1.0 -pname K -nname CL -neutral
fi

