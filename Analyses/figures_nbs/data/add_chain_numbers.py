from pymol import cmd
import sys

in_pdb = sys.argv[1] # Input PDB path
out_pdb = sys.argv[2] # Output PDB path

cmd.load(in_pdb,"MyProtein")
Chains = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H']
for k in range(len(Chains)):
	cmd.alter("chain "+Chains[k],"resi=str(int(resi)-"+str(32*k)+")")

cmd.set("retain_order",1)
cmd.save(out_pdb ,"MyProtein")
