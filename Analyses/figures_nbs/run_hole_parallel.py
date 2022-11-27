import os
import re
import sys
import pickle
import numpy as np
import subprocess
from joblib import Parallel, delayed

def encode_hole_inputs(pdbfile):
    """
    Enconde a sample HOLE input file content and return its string chain
    """
    file_content = ("!HOLE input file",
                    f"coord {pdbfile}",
                    "radius ./hole2/rad/simple.rad",
                    "ignore hoh",
                    "capsule",
                    "endrad 15.")
    
    inputs_encoded = '\n'.join(file_content).encode('utf-8')
    
    return inputs_encoded

def run_hole(inputs_encoded):
    """
    Run hole executable via subprocess with an encoded HOLE input file
    Output file content is returned as an enconded string chain
    """
    command = ['./hole2/exe/hole']
    
    status = 1
    n_trial = 0
    n_max = 4
    while (status != 0 and n_trial <= n_max):
        child = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE)
        outputs_encoded = child.communicate(input=inputs_encoded)[0]
        status = child.returncode
        
        if status != 0:
            n_trial += 1
            mssg = "ERROR. HOLE failed. I will try maximum ", n_max,"times"
            print(mssg, file=sys.stderr)
            
    if status == 0:
        return outputs_encoded
    
    else:
        return 'Failed'.encode('utf-8')
    
    
def process_hole_output(hole_output, key_name = 'pdbname'):
    """
    Compute HOLE conductance estimates and pore lumen dimensions/features.
    If HOLE fails to produce an output, then a null dictionay is returned.
    """
    try:
        for line in hole_output:
            if re.search(r'"Atomic" length of channel',line):
                pore_length = float(line.split()[4]) # [Angstroms]
                
            elif re.search(r'TAG',line):
                # All conductance estimates are in [nano-Siemens]
                x = line.split()
                VDW_Rmin = float(x[3]) # [Angstroms]
                Gmacro = 0.001*float(x[5])
                Gpred_Rmin = 0.001*float(x[8])
                Gpred_Lenght = 0.001*float(x[9])
                Gpred_AvgEPot = 0.001*float(x[10])

                HOLE_dimensions = (VDW_Rmin,pore_length)
                HOLE_conductance_estimates = (Gmacro,Gpred_Rmin,Gpred_Lenght,Gpred_AvgEPot)

        return {key_name:(HOLE_dimensions,HOLE_conductance_estimates)}
    
    except Exception:
        print("ERROR. HOLE failed for:", key_name)
        
        HOLE_dimensions = (0,0)
        HOLE_conductance_estimates = (0,0,0,0)
        
        return {key_name:(HOLE_dimensions,HOLE_conductance_estimates)}
    
def extract_hole_results(pdb):
    inputs_encoded = encode_hole_inputs(pdb)
    outputs_encoded = run_hole(inputs_encoded)
    outputs_decoded = outputs_encoded.decode("utf-8").split('\n')
    
    key_name = os.path.basename(pdb).split('.')[0]
    outputs_dict = process_hole_output(outputs_decoded, key_name=key_name)
    
    return outputs_dict

if __name__ == "__main__":
    from time import process_time
    
    pdbfiles = [x.strip() for x in sys.stdin.readlines()]
    print("INFO. You inputted ", len(pdbfiles), " PDB files to process")
    
    n_cores = int(sys.argv[1])
    output_dir = sys.argv[2]
    
    #tic = process_time()
    results = Parallel(n_jobs=n_cores)(delayed(extract_hole_results)(pdb) for pdb in pdbfiles)
    results_dict = {k: v for x in results for k, v in x.items()}
    
    results_dict.update({'metadata' : {
            'pdb_filename' : ('HOLE_dimensions', 'HOLE_conductance_estimates'),
            'HOLE_dimensions' : ('VDW_Rmin', 'pore_length'),
            'HOLE_conductance_estimates' : ('Gmacro','Gpred_Rmin','Gpred_Lenght','Gpred_AvgEPot')}
                        })
    
    pickle_file = os.path.join(output_dir, "md_100ns.hole.pickle")
    with open(pickle_file, 'wb') as f:
        pickle.dump(results_dict, f)
    
    #toc = process_time()
    #elapsed_time = np.round((toc-tic)/60, decimals=3)
    #print(f"INFO. Total elapsed time [m]: {elapsed_time}")
