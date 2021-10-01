from pymol import cmd
import sys

pdb1 = sys.argv[1] # PDB docking model
pdb2 = sys.argv[2] # PDB md_100ns model
wza_pdb = "data/md_selected_models/wzaD4-consensus-full.pdb"

cmd.load(pdb1,"docking")
cmd.load(pdb2,"md_100ns")
cmd.load(wza_pdb,"wza")

cmd.align("docking", "wza")
cmd.align("md_100ns", "wza")

cmd.hide("spheres","all")
cmd.show("cartoon","all")
cmd.show("sticks",'resname tyr+cys')
cmd.set("orthoscopic",1)
cmd.bg_color('white')

