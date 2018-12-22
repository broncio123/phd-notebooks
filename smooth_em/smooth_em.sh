#! bin/bash

sem_stages=(
	"PR-POPC_SC+BB_Protein_Cbonds" 
	"PR-SC+BB_Protein_Cbonds" 
	"PR-BB_Protein_Cbonds" 
	"No-PR_Cbonds" 
	"No-PR_Ubonds"
	)

init_file="ionise.pdb"
for stage in "${sem_stages[@]}"; do
    gmx_mpi grompp -f mdpf/sem_${stage}.mdp -c $init_file -p topol.top -o em/sem_${stage}.tpr
    mpirun -np 20 mdrun_mpi -s em/sem_${stage}.tpr -deffnm em/sem_${stage}
    wait
    gmx_mpi energy -f em/sem_${stage}.edr -o em/energy_sem_${stage}.xvg
    init_file="em/sem_${stage}.gro"
done