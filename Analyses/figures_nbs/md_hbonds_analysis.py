def get_original_chain(x, peptide_length = 32):
    """
    Return the parent chain number in a homomeric peptide assembly for given a residue number
    """
    N = peptide_length
    if x%N == 0:
        return x//N -1
    else:
        return x//N

def get_original_resnum(x, peptide_length = 32):
    """
    Return the residue number in a homomeric peptide assembly for given a residue number
    """
    N = peptide_length
    if x%N == 0:
        return N
    else:
        return x%N

def find_resid(id, protein_atoms):
    """
    Return original residue number and chain annotations in a homomeric peptide assembly for a given atom id
    """
    f = get_original_chain
    g = get_original_resnum
    return [(str(f(atom.resid)), str(g(atom.resid)), atom.resname) for atom in protein_atoms if id == atom.id]

def intra_to_inter_hbond_counts(counts, u, protein_atoms):
    """
    Compute what percentage of the time intrachain and interchain sidechain H-bonds are found
    """
    interchain_counts = []
    intrachain_counts = []
    for donor, hydrogen, acceptor, count in counts:
        # get biopython objects
        d, h, a = u.atoms[donor], u.atoms[hydrogen], u.atoms[acceptor]
        # get parent labels
        donor_chain, donor_resn, donor_resname = find_resid(d.id, protein_atoms)[0]
        acceptor_chain, acceptor_resn, acceptor_resname = find_resid(a.id, protein_atoms)[0]
        if donor_chain == acceptor_chain:
            intrachain_counts.append(count)
        else:
            interchain_counts.append(count)

    intrachain_frequency = sum(intrachain_counts)/sum(counts.T[-1])
    interchain_frequency = sum(interchain_counts)/sum(counts.T[-1])
    
    return intrachain_frequency, interchain_frequency

def get_leading_hbonds(hbonds_object, protein_atoms):
    """
    Return list of annotated leading H-bonds in a homomeric peptide assembly
    """
    logging.info('Running unrefined H-bonds analysis for protein-protein as donor-acceptor combination')
    hbonds = hbonds_object
    hbonds.hydrogens_sel = hbonds.guess_hydrogens("protein")
    hbonds.acceptors_sel = hbonds.guess_acceptors("protein")
    hbonds.run()

    logging.info('Finding H-bond types with the highest counts')
    hbond_counting_raw = hbonds.count_by_type()
    hbond_counting_sorted = sorted(hbond_counting_raw, key=lambda x: int(x[-1]), reverse=True)
    hbond_counting_sorted
    logging.info('Found %s H-bond types in your trajectory', len(hbond_counting_sorted))
 
    # NOTE: Here I'm taking the top 16, but this is arbitrary
    logging.info('Filtering names of contributing residues')
    data_destilated = hbond_counting_sorted[:16]
    set_donors_acceptors = set([(x[0].split(':')[0], x[1].split(':')[0]) for x in data_destilated])
    set_donors_acceptors
    logging.info('Found a total of %s unique donor-acceptor residue name pairs', len(set_donors_acceptors))

    logging.info('Listing pairs of H-donor and Acceptor residue numbers')
    hdonor_acceptor_list = []
    for data in hbonds.count_by_ids():
        hdonor_atom_id, _, acceptor_atom_id, _ = data
        hdonor = find_resid(hdonor_atom_id, protein_atoms)
        acceptor =  find_resid(acceptor_atom_id, protein_atoms)
        hdonor_acceptor = hdonor + acceptor
        hdonor_acceptor_list.append(('-'.join(hdonor[0]), '-'.join(acceptor[0])))

    logging.info('Splitting list into intrachain and interchain H-bonds') 
    interchain_hbonds = []
    intrachain_hbonds = []
    for hdonor, acceptor in list(set(hdonor_acceptor_list)):
        hdonor_chain, hdonor_resnum, hdonor_resname = hdonor.split('-')
        acceptor_chain, acceptor_resnum, acceptor_resname = acceptor.split('-')
        if hdonor_chain == acceptor_chain:
            intrachain_hbonds.append(('-'.join([hdonor_resnum, hdonor_resname]),
                                      '-'.join([acceptor_resnum, acceptor_resname])))
        else:
            interchain_hbonds.append(('-'.join([hdonor_resnum, hdonor_resname]),
                                      '-'.join([acceptor_resnum, acceptor_resname])))

    logging.info('Filtering leading INTER-chain H-bond residue pairs')
    X = set(interchain_hbonds)
    Y = []
    for pair in X:
        hdonor, acceptor = pair
        pair_clean = (hdonor.split('-')[1], acceptor.split('-')[1])
        if pair_clean in set_donors_acceptors:
            Y.append(pair)
    interchain_hbonds_leads = Y

    logging.info('Filtering leading INTRA-chain H-bond residue pairs') 
    X = set(intrachain_hbonds)
    Y = []
    for pair in X:
        hdonor, acceptor = pair
        pair_clean = (hdonor.split('-')[1], acceptor.split('-')[1])
        if pair_clean in set_donors_acceptors:
            Y.append(pair)
    intrachain_hbonds_leads = Y

    hbonds_leads = set(interchain_hbonds_leads + intrachain_hbonds_leads)
    return hbonds_leads


def get_leading_hbonds_timeseries(hbonds_object, hbonds_leads):
    """
    Extract timeseries for leading inter/intra-chain H-bonds in a homomeric peptide assembly
    """
    N_residues, N_chains = (32, 8)
    extract_resn = lambda x:int(x.split('-')[0])
    
    hbonds = hbonds_object
    HBA_data = {}
    for i in range(len(list(hbonds_leads))):
        hdonor_res, acceptor_res = list(hbonds_leads)[i]
        logging.info('Running analysis for %s', list(hbonds_leads)[i])
        
        resid1, resid2 = (extract_resn(hdonor_res), extract_resn(acceptor_res))
        donor_selection = 'not backbone and (resid '+ \
                       ' '.join(map(str, N_residues*np.arange(N_chains) + resid1))+')'
        acceptor_selection = 'not backbone and (resid '+ \
                         ' '.join(map(str, N_residues*np.arange(N_chains) + resid2))+')'
        
        X = hbonds.guess_hydrogens(donor_selection)
        X = X.split(' or ')
        Y = [x.split(' and ')[-1].strip(')') for x in X]
        hydrogen_names = ' or '.join(Y) if len(Y) > 1 else Y[0]
        
        X = hbonds.guess_acceptors(acceptor_selection)
        X = X.split(' or ')
        Y = [x.split(' and ')[-1].strip(')') for x in X]
        acceptor_names = ' or '.join(Y) if len(Y) > 1 else Y[0]
        
        hydrogens_sel = donor_selection+' and ('+hydrogen_names+')'
        acceptors_sel =  acceptor_selection+' and ('+acceptor_names+')'
        logging.info('Selected DONORS: %s', hydrogens_sel)
        logging.info('Selected ACCEPTORS: %s', acceptors_sel)
        
        hbonds.hydrogens_sel = hydrogens_sel
        hbonds.acceptors_sel = acceptors_sel
        tic = process_time()
        hbonds.run()
        toc = process_time()
        logging.info('Elapsed time (secs): %s', (toc-tic))
        
        HBA_data[list(hbonds_leads)[i]] = {'ids': hbonds.count_by_ids(),
                                           'type': hbonds.count_by_type(),
                                           'time': hbonds.count_by_time()}
    
    return HBA_data

if __name__ == "__main__":
    import os
    import sys
    import pickle
    import coloredlogs, logging
    import argparse
    import numpy as np
    from time import process_time
    import MDAnalysis as mda
    from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis
    
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
                        help="Set H-bond analysis output filename")
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
        logging.basicConfig(filename=os.path.join(workdir,'HBA.log'),
                            filemode='w', 
                            format='%(name)s - %(levelname)s - %(message)s')

    # Load trajectory file along with topology
    tic = process_time() 
    u = mda.Universe(topology_filepath, trajectory_filepath, in_memory=True)
    toc = process_time()
    logging.info('Loaded your trajectory and topology.')
    logging.info('Elapsed time (secs): %s', (toc-tic))

    # Select protein atoms
    protein_atoms = list(u.select_atoms("protein"))
    logging.info('There are %s atoms in your protein selection.', str(len(protein_atoms)))
    
    # Set HBA object
    logging.info('Initialising HydrogenBondAnalysis object instance with your universe')
    hbonds = HydrogenBondAnalysis(universe=u)
    
    logging.info('Working out Leading Hydrogen-bonds in your trajectory.')
    tic_total = process_time() 
    leading_hbonds = get_leading_hbonds(hbonds, protein_atoms)
    toc_total = process_time() 
    logging.info('Found %s Leading Hydrogen-bonds.', str(len(leading_hbonds)))
    logging.info('Elapsed time (mins): %s', (toc_total-tic_total)/60)
    
    logging.info('Extracting H-bond timeseries for Leading Hydrogen-bonds.')
    tic_total = process_time() 
    timeseries = get_leading_hbonds_timeseries(hbonds, leading_hbonds)
    toc_total = process_time() 
    logging.info('Elapsed TOTAL time (mins): %s', (toc_total-tic_total)/60)
    
    #Save the data
    output_filepath = os.path.join(workdir, output_filename)
    with open(output_filepath,'wb') as f:
        pickle.dump(timeseries, f)
    logging.info('Saved your H-bond Analysis data in %s', output_filepath)
