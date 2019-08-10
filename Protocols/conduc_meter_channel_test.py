#!/usr/bin/env python

# Note: Modified version 18 Nov 2016
import sys, os, numpy, itertools, subprocess
from conduc_meter_library import *

fname   = sys.argv[1] # Generic name of simulation files
option  = sys.argv[2] # Output data option (Print or Save)

info    = subprocess.Popen(['tail','-1',fname+'.gro'],stdout=subprocess.PIPE)
box_dim = map(float,info.stdout.read().split()) # Box dimensions [nm ,nm, nm]
Lz      = box_dim[-1]   # Box Z-length [nm]

# Extract atom indexes and coordinates
indexes_and_coordinates(fname)
# End coordinates of bilayer leaflets [nm]
Z_min, Z_max    = find_P8end_coords(fname)
# Determine relative permeations
permeations = instant_charge(fname, Lz, Z_min, Z_max)

N       = len(permeations['K'])
if option == 'save':
        outfile = open(fname+'.cm_flux', 'w')
        for n in range(N):
                time       = permeations['K'][n][0]
                icharge_K  = permeations['K'][n][1]
                icharge_CL = permeations['CL'][n][1]
                outfile.write( str(time)+'\t'+str(icharge_K)+'\t'+str(icharge_CL)+'\n' )
        outfile.close()
elif option == 'print':
        for n in range(N):
                time       = permeations['K'][n][0]
                icharge_K  = permeations['K'][n][1]
                icharge_CL = permeations['CL'][n][1]
                print time, icharge_K, icharge_CL
else:
        print "Provide a valid option: 'save' or 'print'"