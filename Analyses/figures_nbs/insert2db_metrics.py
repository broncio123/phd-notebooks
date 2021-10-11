import sys
import os
import numpy
import operator
import subprocess
import json
import concurrent.futures
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

# Import modules from folder
modules_path = "/home/ba13026/mpmodeling/analysis/"
if modules_path not in sys.path:
    sys.path.append(modules_path)

import analyse_HOLE
from setup_db_metrics import Base, Tags, Pore_Dimensions, Radii_of_Gyration

def get_tags(model_pdb):
    idents = operator.itemgetter(*[0,1,2,-1])(model_pdb.split('/')[9:])
    mutant, group, model_name, frame = idents
    frame = frame[:-4]
    return mutant, group, model_name, frame
     
def get_Rg_components(model_pdb):
    try:
        u = mda.Universe(model_pdb)
        protein = u.select_atoms('protein and backbone')
        protein_mass = numpy.sum(protein.masses)
        protein_inertia = protein.moment_of_inertia() # tensor of inertia
        f = lambda x : numpy.sqrt(x/float(protein_mass))
        protein_Rg_n = [f(protein_inertia[i][i]) for i in range(3)]
        Rg_x, Rg_y, Rg_z = protein_Rg_n
        return Rg_x, Rg_y, Rg_z
    except:
        return 'Fail'

def get_HOLE_Rmin(model_pdb):
    dir_path = os.getcwd()
    tmp_wd = os.path.dirname(model_pdb)
    os.chdir(tmp_wd)
    try:
        HOLE_dimensions,HOLE_conductance_estimates = analyse_HOLE.hole(os.path.basename(model_pdb))
        os.chdir(dir_path)
        return HOLE_dimensions[0]
    except:
        os.chdir(dir_path)
        return 'Fail'
        
def get_channel_length(model_pdb):
    try:
        model_ampal = isambard_dev.ampal.convert_pdb_to_ampal(model_pdb)
        get_chain_Zcoords = lambda chain : [atom.z for atom in chain.get_atoms()]
        get_length = lambda chain : max(get_chain_Zcoords(chain)) - min(get_chain_Zcoords(chain))
        chains_lengths = list(map(get_length, model_ampal))
        return numpy.mean(chains_lengths)
    except:
        return 'Fail'

def process_model(n):
    model_pdb = param_list[n]
    #####################################
    # Model identifiers
    mutant, group, model_name, frame = get_tags(model_pdb)
    model_tags  = Tags(
                        mutant = mutant, 
                        group = group, 
                        pdb_name = model_name,
                        frame = frame
                        )
    session.add(model_tags)    
    #####################################
    # Radius of gyration decomposition
    data = get_Rg_components(model_pdb)
    model_Rgs = Radii_of_Gyration(
        Rg_x = data[0],
        Rg_y = data[1],
        Rg_z = data[2],
        tag  = model_tags
        )
    session.add(model_Rgs)
    #####################################
#     # Pore dimensions
    pore_Rmin = get_HOLE_Rmin(model_pdb)
    pore_length = get_channel_length(model_pdb)
    model_pore_dimensions = Pore_Dimensions(
        pore_Rmin = pore_Rmin,
        pore_length = pore_length,
        tag = model_tags
        )
    session.add(model_pore_dimensions)
    #####################################
    # COMMIT CHANGES TO DATABASE  
    session.commit()

#####################################
# Parallel Process Execution
def main():
    model_n = list(range(len(param_list)))
    with concurrent.futures.ProcessPoolExecutor(max_workers = ncores) as executor:
        executor.map(process_model, model_n)

if __name__ == '__main__':
    dbfile  = sys.argv[1] # Database filename
    param_json = sys.argv[2] # Dictionary with mutant structures info
    ncores = int(sys.argv[3]) # Number of cores

    # Extract info from dictionary
    with open(param_json, 'r') as fp:
        param_list = json.load(fp)
    fp.close()
    
    # Create engine and bind it to current session
    engine = create_engine('sqlite:///'+dbfile)
    Base.metadata.bind = engine
    DBSession = sessionmaker(bind=engine)
    session = DBSession()

    main()

