#!/usr/bin/env python

import sys, os, numpy, itertools, subprocess

def indexes_and_coordinates(fname):
        # Test if index files already exist
        Atoms = ['P8', 'K', 'CL'] # Define atom names
        # Make individual index and trajectories files
        for atom in Atoms:
                print ('TESTING IF INDEX FILE FOR ATOM %s IS PRESENT IN CURRENT DIRECTORY ...\n'% (atom))
                if os.path.isfile(atom+'-Indexes_'+fname+'.ndx') == False:
                        print ('Index file does not exist for atom: '+atom+'\n')
                        print ('Index file will be created using GROMACS tool: make_ndx ...\n')
                        gmx_code0       = ['gmx','make_ndx','-f',fname+'.gro','-o',atom+'-Indexes_'+fname+'.ndx']
                        proc            = subprocess.Popen(gmx_code0, stdin=subprocess.PIPE)
                        proc.communicate(input = 'del 0-50\n a '+atom+'\n q\n')
                        proc.wait()
                else:
                        print('INDEX FILE %s ALREADY EXISTS! ...\n'% (atom+'-Indexes_'+fname+'.ndx'))

                print ('TESTING Z-COORDINATE FILE FOR ATOM %s IS PRESENT IN CURRENT DIRECTORY ...\n'% (atom))
                if os.path.isfile(atom+'-Zcoord_'+'PBC_'+fname+'.xvg') == False:
                        print ('File with Z-coordinate does not exist for atom: '+atom+'\n')
                        print ('Z-coordinates will be extracted using GROMACS tool: traj ...\n')
                        gmx_code1       = ['gmx','traj','-f',fname+'.xtc','-s',fname+'.tpr','-n',atom+'-Indexes_'+fname+'.ndx','-nox','-noy','-ox',atom+'-Zcoord_'+'PBC_'+fname+'.xvg']
                        subprocess.call(gmx_code1)
                        print ('Coordinates successfully extracted :)\n')
                        print ('Keep in mind that these coordinates are wrapped (not corrected by PBCs)\n')
                else:
                        print('Z-COORDINATE FILE %s ALREADY EXISTS!\n'% (atom+'-Zcoord_'+'PBC_'+fname+'.xvg') )

# Function to correct displacement by PBC
def pbc_distance(x,y,Lz):
    d = y - x
    if abs(d) <= Lz/float(2):
        return d
    else:
        return numpy.sign(d)*(abs(d)-Lz)

def find_P8end_coords(fname):
        #Average position of P8 atoms
        infile  = open('P8-Zcoord_'+'PBC_'+fname+'.xvg','r')    # Open file
        lines   = infile.readlines()    # Read lines
        N       = len(lines)    # Number of lines
        infile.close()
        mZ_P8_UP        = []
        mZ_P8_DOWN      = []
        for n in range(N):
                if  ("@" not in lines[n].rstrip()) and ("#" not in lines[n].rstrip()):
                        Z_P8    = map(float, lines[n].rstrip().split()[1:]) # Z-ccordinates P8 atoms
                        Z0      = numpy.mean(Z_P8)
                        mZ_P8_UP.append( numpy.mean([z for z in Z_P8 if z > Z0]) )
                        mZ_P8_DOWN.append( numpy.mean([z for z in Z_P8 if z < Z0]) )
        P8_UP, P8_DOWN = numpy.mean(mZ_P8_UP), numpy.mean(mZ_P8_DOWN)
        return P8_DOWN, P8_UP # P8-end coordinates [nm]

def instant_charge(fname, Lz, Z_min, Z_max):
        """This code computes both the instant normalised-charge [C/1.6E-19] permeated in time [ps], for each ion species."""
        permeations = {} # Relative permeations
        bthickness  = numpy.abs(Z_max - Z_min)
        for ion in ['K', 'CL']:
                # Load ion indexes
                ndxfile = open(ion+'-Indexes_'+fname+'.ndx')
                ion_ndx = list(itertools.chain.from_iterable([x.rstrip().split() for x in ndxfile.readlines()[1:]]))
                Nc      = len(ion_ndx)  # Number of ions
                ndxfile.close()

                # Trajectory file
                infile  = open(ion+'-Zcoord_'+'PBC_'+fname+'.xvg', 'r')
                lines   = infile.readlines()    # Read lines
                N       = len(lines)    # Number of lines
                infile.close()

                Time    =       []      # Time [ps]
                ICharge =       []      # 'Relative' permeation events (Instant charge)

                for n in range(N-1):
                        if  ("@" not in lines[n].rstrip()) and ("#" not in lines[n].rstrip()):
                                Time.append( float(lines[n+1].rstrip().split()[0]) )
                                # Z-coordinates of ions at time 't' (Current frame)
                                Z_t1      = map(float, lines[n].rstrip().split())
                                # Z-coordinates of ions at time 't+dt' (Next frame)
                                Z_t2      = map(float, lines[n+1].rstrip().split())
                                # List all pairs of coordinates per ion at times 't' and 't+dt'
                                Pairs   = [x for x in itertools.izip(Z_t1, Z_t2)]
                                # Compute Forward Displacement, correcting by PBCs
                                dZ      = [pbc_distance(x,y,Lz) for x,y in Pairs]
                                # Check if ions inside channel volume at time 't+dt'
                                State   = ['In'  if  Z_min <= x <= Z_max else 'Out' for x in Z_t2]
                                # Filtre displacement values of ions inside channel at time 't+dt'
                                dZ_in   = [dZ[k] if State[k]=='In' else 0 for k in range(len(dZ))]
                                # Compute 'relative' permeation events
                                ICharge.append( sum(dZ_in)/float(bthickness) ) # Instant Charge [C/1.6E-19]

                permeations[ion] = numpy.array([Time, ICharge]).T
        return permeations