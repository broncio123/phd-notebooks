import numpy as np
import pickle
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis.analysis import distances
from itertools import combinations

mutant = "cWza-Y373C"

# Load trajectory
workdir = "data/cWza-Y373C_conformation1_0337/"
u = mda.Universe(workdir+'md_100ns.tpr', workdir+'md_100ns.xtc',in_memory=True)

# work out residue pair combinations
residue_combinations = list(combinations(range(1,257), 2))

# Define function to search for contacts and output their timeseries
def contacts_within_cutoff(u, group_a, group_b, radius=4.5):
    timeseries = []
    for ts in u.trajectory:
        # calculate distances between group_a and group_b
        dist = contacts.distance_array(group_a.positions, group_b.positions)
        # determine which distances <= radius
        n_contacts = contacts.contact_matrix(dist, radius).sum()
        timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries).T[1]

contacts_data = {}
for res_pair in residue_combinations:
    res_a, res_b = res_pair
    # Define groups - select heavy atoms only, exclude Hydrogens
    group_a = u.select_atoms("resid "+str(res_a)+" and not(name H* or backbone)")
    group_b = u.select_atoms("resid "+str(res_b)+" and not(name H* or backbone)")
    timeseries = contacts_within_cutoff(u, group_a, group_b, radius=4.5)
    if timeseries.sum() > 0:
        contacts_data[res_pair] = timeseries

# Save data
with open("data/contacts_data_"+mutant+".pkl", "wb") as fp:
    pickle.dump(contacts_data, fp)
