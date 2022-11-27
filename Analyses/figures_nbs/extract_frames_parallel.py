import os
import sys
import numpy as np
import mdtraj as md
from time import process_time
from joblib import Parallel, delayed

# Example:
# python extract_frames_parallel.py data/cWza-Y373C_conformation1_0209/ 8

dirname = sys.argv[1]
n_cores = int(sys.argv[2])

tic = process_time()
traj_file = os.path.join(dirname, 'md_100ns.xtc')
top_file = os.path.join(dirname, 'md_100ns.gro')
trajectory = md.load(traj_file, top=top_file)
protein_indices = trajectory.topology.select("protein")

framesdir = os.path.join(dirname, 'frames')
if os.path.isdir(framesdir):
    pass
    print("INFO: Dir 'frames' exists.")
else:
    os.mkdir(framesdir)
    print("INFO: Created 'frames' dir.")

def extract_protein_frame(i, trajectory, protein_indices, outdir):
    frame = trajectory[i]
    protein_frame = frame.restrict_atoms(atom_indices=protein_indices)
    pdb_output = os.path.join(outdir, 'frame'+f'{i:04d}'+'.pdb')
    protein_frame.save_pdb(filename=pdb_output)
    
print(f"INFO: Ypu have {trajectory.n_frames} frames in your trajectory: {dirname}")
f = lambda i:extract_protein_frame(i, trajectory, protein_indices, framesdir)
p = Parallel(n_jobs=n_cores)(delayed(f)(i) for i in range(trajectory.n_frames))

toc = process_time()
elapsed_time = np.round((toc-tic)/60, decimals=3)
print(f"INFO: Total elapsed time [m]: {elapsed_time}")
