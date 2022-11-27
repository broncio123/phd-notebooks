import sys
import os
import coloredlogs, logging
import argparse
import numpy as np
from time import process_time
import MDAnalysis as mda
from MDAnalysis.analysis import hole2

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
                    help="Set summary file from HOLE conductance analysis output filename")
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
    logging.basicConfig(filename=os.path.join(workdir,'hole.log'),
                        filemode='w', 
                        format='%(name)s - %(levelname)s - %(message)s')

# Load trajectory file along with topology
tic = process_time() 
u = mda.Universe(topology_filepath, trajectory_filepath, in_memory=True)
toc = process_time()
logging.info('Loaded your trajectory and topology.')
logging.info('Elapsed time (secs): %s', (toc-tic))

# HOLE settings
hole_path = '/data/dragon000/sanjuan/bristol/cwza/hole2/exe/hole' # absolute path to executable

outdir = os.path.join(workdir, output_filename)
if os.path.isdir(outdir):
    logging.warning("Dir already exists: %s", outdir)
else:
    logging.info("I'll create a new dir for: %s", outdir)
    os.mkdir(outdir)

tic_total = process_time()
hole_analysis = hole2.HoleAnalysis(u, executable=hole_path, prefix = outdir+'/')
hole_analysis.run()
toc_total = process_time()
logging.info('Finished HOLE conductance analysis for model: %s', workdir)
logging.info('Elapsed TOTAL time (mins): %s', (toc-tic)/60)
