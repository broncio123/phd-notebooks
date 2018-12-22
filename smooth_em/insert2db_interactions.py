import sys
import numpy
import operator
import subprocess
import json
import concurrent.futures
import isambard_dev
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from setup_db_interactions import  Json, Pdb, Interhelix_Interactions, Chain2Complex_BB_distance, HOLE_Output, Base
import analyse_HOLE


def interaction_direction(chain_combination):
    chainc_ccwise = ['AB','BC','CD','DE','EF','FG','GH','HA']
    chainc_cwise  = ['AH','HG','GF','FE','ED','DC','CB','BA']
    chainc_intrachn = ['AA','BB','CC','DD','EE','FF','GG','HH']
    if chain_combination in chainc_cwise:
        return 1
    elif chain_combination in chainc_ccwise:
        return -1
    elif chain_combination in chainc_intrachn:
        return 0

def get_OH_atoms(file):
	"""Get all OH-atoms per H-bond in PDB, in residue-number/OH-atom format plus chain direction"""
	try:
		p = isambard_dev.ampal.convert_pdb_to_ampal(file)
		hbonds = isambard_dev.interactions.find_hydrogen_bonds(p)
		# Find all H-bonds and select those between sidechain atoms
		sc_hbonds_raw = [hb for hb in hbonds if hb.is_sidechain_sidechain == True]
		sc_hbonds_reduced = []
		for hb in sc_hbonds_raw:
			donor_H = [hb.donor.ampal_parent.id , hb.donor.res_label]
			acceptor_O = [hb.acceptor.ampal_parent.id , hb.acceptor.res_label]
			direction_HO = hb.donor.unique_id[0]+hb.acceptor.unique_id[0]
			hb_reduced = donor_H+acceptor_O+[direction_HO]
			data = json.dumps(hb_reduced)
			sc_hbonds_reduced.append( data )
	except:
		sc_hbonds_reduced = 'NoFile'
	return sc_hbonds_reduced

def get_KIHs(file):
    """Get all KIHs in PDB, in residue number format plus chain direction"""
    try:
        p = isambard_dev.ampal.convert_pdb_to_ampal(file)
        kihs_raw = isambard_dev.add_ons.knobs_into_holes.find_kihs(p)
        kihs_reduced = []
        for kih in kihs_raw:
            knob_data = ''.join(kih.knob.unique_id).rstrip()
            hole_data = [''.join(kih.hole[x].unique_id).rstrip() for x in range(len(kih.hole))]
            kih_direction = knob_data[0]+hole_data[0][0]
            kih_reduced = [knob_data[1:]]+[s[1:] for s in hole_data]+[kih_direction]
            data = json.dumps(kih_reduced)
            kihs_reduced.append(data)
    except:
        kihs_reduced = 'NoFile'
    return kihs_reduced

def get_HOLE(file):
	try:
		HOLE_dimensions,HOLE_conductance_estimates = analyse_HOLE.hole(file)
		data = HOLE_dimensions+HOLE_conductance_estimates
		return json.dumps(data)
	except:
		return 'NoFile'

def get_COM_bb_distance(file):
	try:
		p = isambard_dev.ampal.convert_pdb_to_ampal(file)
		ccom = p.centre_of_mass
		n_chains = len(p.sequences)
		com_dd_distances = []
		for k in range(n_chains):
			com_dd_distances.append(numpy.linalg.norm(p[k].centre_of_mass - ccom))
		return json.dumps(com_dd_distances)
	except:
		return 'NoFile'

def data_from_trajectory(path, n_frames, get_data):
	data_per_frame = []
	for n in range(n_frames+1):
		file = path+'frame_'+str(n)+'.pdb'
		data_per_frame.append(get_data(file))
	return data_per_frame


if __name__ == '__main__':
	dbfile	= sys.argv[1] # Database filename
	pdblist = sys.argv[2]  # List of PDB files
	ncores = int(sys.argv[3]) # Number of cores

	# List of paths to pfb frames
	pdb_paths_original = [x.rstrip()+'/complex/mdf/prmd' for x in open(pdblist,'r').readlines()]
	# Extract PDB names
	pdbnames = [x.rstrip().split('/')[2] for x in open(pdblist,'r').readlines()]
	# Create soft links to frame directories
	for k in range(len(pdb_paths_original)):
		cmd = ['ln', '-s', pdb_paths_original[k], './prmd_'+pdbnames[k]]
		subprocess.call(cmd)
	# List of modified paths after creating soft links to folders
	pdb_paths = ['./prmd_'+pdbnames[k]+'/' for k in range(len(pdbnames))]

	# Create engine and bind it to current session
	engine = create_engine('sqlite:///'+dbfile)
	Base.metadata.bind = engine
	DBSession = sessionmaker(bind=engine)
	session = DBSession()

	def process_model(n):
		## Define PDB name	
		model  = Pdb(pdb_name = pdbnames[n])
		session.add(model) 

		##  Inter-helix interactions per frame in MD trajectory: H-bonds and KIHs
		n_frames = 100
		path = pdb_paths[n]
		hbonds_trajectory = data_from_trajectory(path,n_frames,get_OH_atoms)
		kihs_trajectory	= data_from_trajectory(path,n_frames,get_KIHs)
		model_interactions = Interhelix_Interactions(hbonds=hbonds_trajectory,kihs=kihs_trajectory,pdb=model)
		session.add(model_interactions)

		## HOLE dimensions [Angstroms] and conductance estimates [nS]
		hole_trajectory = data_from_trajectory(path,n_frames,get_HOLE)
		model_HOLE = HOLE_Output(HOLE_data=hole_trajectory,pdb=model)	
		session.add(model_HOLE)
		
		## Chain-BB-COM to Assembly-BB-COM distances per frame in MD trajectory
		bb_distances_trajectory = data_from_trajectory(path,n_frames,get_COM_bb_distance)
		model_COM_distances = Chain2Complex_BB_distance(bb_distances=bb_distances_trajectory,pdb=model)
		session.add(model_COM_distances)

		# COMMIT CHANGES TO DATABASE  
		session.commit()

	model_n = list(range(len(pdb_paths)))
	def main():
		with concurrent.futures.ProcessPoolExecutor(max_workers=ncores) as executor:
			executor.map(process_model,model_n)

	main()

