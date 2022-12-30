def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries).T[1]
    
if __name__ == "__main__":
    import os
    import sys
    import pickle
    import coloredlogs, logging
    import argparse
    import numpy as np
    from time import process_time
    import MDAnalysis as mda
    from MDAnalysis.analysis import contacts
    from MDAnalysis.analysis import distances
    from itertools import combinations
    
    # Initilise argparser
    parser = argparse.ArgumentParser(description='This is a demo.')
    parser.add_argument("-p",
                        "--topology",
                        dest="topology_filepath",
                        help="Set topology filepath")
    parser.add_argument("-x",
                        "--trajectory",
                        dest="trajectory_filepath",
                        help="Set trajectory filepath")
    parser.add_argument("-o",
                        "--output",
                        dest="output_filename",
                        help="Set noncovalent contacts analysis output filename")
    parser.add_argument("-v",
                        "--verbose",
                        help="Increase output verbosity",
                        action="store_true")
    
    # Input/Output files
    args = parser.parse_args()
    
    topology_filepath = args.topology_filepath
    trajectory_filepath = args.trajectory_filepath
    output_filename = args.output_filename
    
    # Configure logger
    workdir = os.path.dirname(trajectory_filepath)
    if args.verbose:
        coloredlogs.install(level='DEBUG')
        logging.basicConfig(level=logging.DEBUG)
        logging.basicConfig(filename=os.path.join(workdir,'contacts.log'),
                            filemode='w', 
                            format='%(name)s - %(levelname)s - %(message)s')

    # Load trajectory file along with topology
    tic = process_time() 
    u = mda.Universe(topology_filepath, trajectory_filepath, in_memory=True)
    toc = process_time()
    logging.info('Loaded your trajectory and topology.')
    logging.info('Elapsed time (secs): %s', (toc-tic))

    # Work out residue pair combinations
    n_residues = 32
    n_chains = 8
    residue_combinations = list(combinations(range(1,n_residues*n_chains+1), 2))

    # Loop over all residue combination pairs
    logging.info('Extracting Non-covalent contacts for your input trajectory and topology. Default cut-off = 4.5 (Å).')
    tic_total = process_time() 
    
    contacts_data = {}
    for res_pair in residue_combinations:
        res_a, res_b = res_pair
        # Define groups - select Side-Chain heavy atoms only, exclude Hydrogens
        group_a = u.select_atoms("resid "+str(res_a)+" and not(name H* or backbone)")
        group_b = u.select_atoms("resid "+str(res_b)+" and not(name H* or backbone)")
        timeseries = contacts_within_cutoff(u, group_a, group_b)
        # Ignore null contacts
        if timeseries.sum() > 0:
            contacts_data[res_pair] = timeseries
            
    toc_total = process_time()
    logging.info('Elapsed TOTAL time (mins): %s', (toc_total-tic_total)/60)
    
    #Save the data
    output_filepath = os.path.join(workdir, output_filename)
    with open(output_filepath,'wb') as f:
        pickle.dump(contacts_data, f)
    logging.info('Saved your Non-covalent Contacts Analysis data in %s', output_filepath)