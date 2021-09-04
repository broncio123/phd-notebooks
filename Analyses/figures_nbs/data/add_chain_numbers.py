from pymol import cmd
import sys

pdb_in = sys.argv[1] # Input PDB path
pdb_out = sys.argv[2] # Output PDB path

cmd.load(pdb_in,"MyProtein")
N_atoms = cmd.count_atoms("MyProtein")

Chains = ['A', 'B','C','D','E','F','G','H']

atoms_per_chain = N_atoms/len(Chains)

for k in range(len(Chains)):
    atom_number_intial = int( 1 + k*atoms_per_chain )
    atom_number_final  = int( (k + 1)*atoms_per_chain )
    selection = "id "+ str(atom_number_intial) + ":" + str(atom_number_final)
    expression = "chain='"+Chains[k]+"'"
    cmd.alter(selection, expression)

cmd.set("retain_order",1)
cmd.save(pdb_out ,"MyProtein")
