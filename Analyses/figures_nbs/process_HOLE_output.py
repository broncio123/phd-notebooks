#!/usr/bin/env python
import sys,subprocess,re, os
from ast import literal_eval

def read_hole_output(filepath):
    """
    Compute HOLE conductance estimates and pore lumen dimensions/features
    NOTE: Had to wrap HOLE code in bash code 'run_hole' to automatically generate HOLE input file
    """
    try:
        hole_lines = open(filepath,'r').readlines()
    except FileNotFoundError:
        print("ERROR: File not found")
        HOLE_dimensions = (None,None)
        HOLE_conductance_estimates = (None,None,None,None)
        
        return(HOLE_dimensions,HOLE_conductance_estimates)

    # Filter HOLE output file
    for l in hole_lines:
        if re.search(r'"Atomic" length of channel',l):
            pore_length = float(l.split()[4]) # [Angstroms]
        elif re.search(r'TAG',l):
            # All conductance estimates are in [nano-Siemens]
            x = l.split()
            VDW_Rmin = float(x[3]) # [Angstroms]
            Gmacro = 0.001*float(x[5])
            Gpred_Rmin = 0.001*float(x[8])
            Gpred_Lenght = 0.001*float(x[9])
            Gpred_AvgEPot = 0.001*float(x[10])
    
    HOLE_dimensions = (VDW_Rmin,pore_length)
    HOLE_conductance_estimates = (Gmacro,Gpred_Rmin,Gpred_Lenght,Gpred_AvgEPot)

    return(HOLE_dimensions,HOLE_conductance_estimates)

if __name__ == "__main__":
    inputfile = sys.argv[1] #HOLE output filepath
    
    HOLE_dimensions,HOLE_conductance_estimates = read_hole_output(inputfile)
    VDW_Rmin,pore_length = HOLE_dimensions
    
    print("Filename: ", inputfile)
    print('HOLE VDW_Rmin [Angstroms]: ', VDW_Rmin)
    print('HOLE Pore_Length [Angstroms]: ', pore_length)
    Gmacro,Gpred_Rmin,Gpred_Length,Gpred_AvgEPot = HOLE_conductance_estimates
    print('HOLE Gmacro [nS]: ', Gmacro)
    print('HOLE Gpred by Rmin [nS]: ', Gpred_Rmin)
    print('HOLE Gpred by Length [nS]: ', Gpred_Length)
    print('HOLE Gpred by Avg EPot [nS]: ', Gpred_AvgEPot)
