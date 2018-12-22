import sys
import numpy
import operator
import json
import itertools
import concurrent.futures
import isambard_dev
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from setup_db_interactions import  Json, Pdb, Interhelix_Interactions, Chain2Complex_BB_distance, HOLE_Output, Base
import analyse_HOLE
############################################################
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
############################################################
def get_data(var, session):
	if var == 'bb_distances':
		return [model.bb_distances for model in session.query(Chain2Complex_BB_distance).all()]
	elif var == 'HOLE_data':
		return [model.HOLE_data for model in session.query(HOLE_Output).all()]
	elif var == 'hbonds':
		return [model.hbonds for model in session.query(Interhelix_Interactions).all()]
	elif var == 'kihs':
		return [model.kihs for model in session.query(Interhelix_Interactions).all()]
	else:
		print("Variable not found. Check database variable defitions.")

##############################################################
def get_time_average(model_data, S_x_reduced_ordered, n_chains):
	n_frames = len(model_data)
	prob_per_traj = []
	for k2 in range(n_frames):
	    x_frame = model_data[k2]
	    prob_per_frame = []
	    for S_x in S_x_reduced_ordered:
	        counter_x = 0
	        for x in x_frame:
	            x_data_raw = json.loads(x)
	            x_direction = [interaction_direction(x_data_raw[-1])]
	            x_atoms = x_data_raw[:-1]
	            x_reduced = json.dumps(x_atoms+x_direction)
	            if x_reduced == S_x:
	                counter_x += 1
	        prob_per_frame.append(counter_x/float(n_chains))
	    prob_per_traj.append(prob_per_frame)
	prob_time_mean =  list(numpy.mean(numpy.array(prob_per_traj).T,axis=1))
	return prob_time_mean

##############################################################
def get_group_average(data, S_x_reduced_ordered, n_chains):
	n_models = len(data)
	prob = []
	for k1 in range(n_models):
	    if 'NoFile' not in data[k1]:
	        model_data = data[k1]
	        prob_time_mean = get_time_average(model_data, S_x_reduced_ordered, n_chains)
	        prob.append(prob_time_mean)
	stats = list(numpy.mean(numpy.array(prob).T,axis=1))
	return stats
#############################################################
def get_S_x(data):
	S_x = set()
	n_models = len(data)
	for k1 in range(n_models):
	    if 'NoFile' not in data[k1]:
	        model_data = data[k1]
	        n_frames = len(model_data)
	        S_x_model = set(itertools.chain.from_iterable(data[k1]))
	        S_x = S_x.union(S_x_model)
	return S_x
##############################################################
def get_S_x_reduced(S_x):
	S_x_reduced = []
	for x_raw in list(S_x):
	    x_data = json.loads(x_raw)
	    x_direction = [interaction_direction(x_data[-1])]
	    x_atoms = x_data[:-1]
	    x_reduced = json.dumps(x_atoms+x_direction)
	    S_x_reduced.append(x_reduced)
	    
	S_x_reduced =  list(set(S_x_reduced))
	return S_x_reduced
#############################################################################################################################
