#!/bin/bash

# Note: 21 Nov 2016
#1. Need to make code for 'create_inp'. Make sure pathways are correct. (Done! 22-11-16)
#2. Need to make code for 'G_stats.py'. Print out averages and std deviations.
#	Select most relevant quantities we want to keep track of.
#3. Test code to death! 

#Note: 22 Nov 2016
#1. Need to rename frame files to feature time-points[ps]
#2. Re-think of how many essays we need to get the std deviation on HOLE calculations
# Keep in in mind that the more essays, the longer the analysis will take.  

function create_inp {
# This function creates ai simple HOLE input file
fname=$1
echo "!HOLE input file
coord ${fname}.pdb
radius ./hole2/rad/simple.rad
ignore hoh
capsule
endrad 15." > ${fname}.inp
}

pdbfile=$1 # PDB file name
# Create a HOLE input file per frame file
create_inp ${pdbfile%.pdb}
# Run HOLE for each frame file (5 independent trials)
./hole2/exe/hole < ${pdbfile%.pdb}.inp > ${pdbfile%.pdb}.txt
# Remove input file
rm -rf ${pdbfile%.pdb}.inp
