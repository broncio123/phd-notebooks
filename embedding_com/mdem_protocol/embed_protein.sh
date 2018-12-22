#!/bin/sh

module load apps/gromacs-5.0.6

outfolder=$1

# Run embedding
gmx_mpi grompp -f ${outfolder}/complex/mdpf/embed.mdp -c  ${outfolder}/complex/for_embedding.gro -p ${outfolder}/complex/for_embedding.top -o ${outfolder}/complex/for_embedding.tpr -n ${outfolder}/complex/for_embedding.ndx

printf "1\n24\n"| mdrun_mpi -s ${outfolder}/complex/for_embedding.tpr -membed ${outfolder}/complex/mdpf/membed.dat -mn ${outfolder}/complex/for_embedding.ndx -mp ${outfolder}/complex/for_embedding.top -deffnm ${outfolder}/complex/confout
# Remove extra generated files
rm ./step*

