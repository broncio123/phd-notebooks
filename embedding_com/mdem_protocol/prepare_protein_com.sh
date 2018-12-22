
module add apps/gromacs-5.0.6
module add boost/gcc-6.1.0/1.64.0

protein_pdb=$1 # Protein input name
protein=$(basename ${protein_pdb%.pdb})
outfolder=$2

mkdir ${outfolder}/bilayer
mkdir ${outfolder}/complex
mkdir ${outfolder}/complex/jobf

DIR=$(dirname $0)
cp -r ${DIR}/templates/embedding/* ${outfolder}/complex/
cp ${DIR}/bilayer_drier.py ${outfolder}/complex/
cp $protein_pdb $outfolder

# Create GROMACS coordinate and topology files for protein
        ## Options: OPLS-AA force field; No ionisation of termini; Add hydrogens
HIS_option=$(bash ~/mpmodeling/find_histype.sh ${outfolder}/${protein}.pdb)
for n in `seq 1 8`; do printf "%s\n" $HIS_option; done | gmx_mpi pdb2gmx -f ${outfolder}/${protein}.pdb  -water none -ff oplsaa -his -o ${outfolder}/confout_pep.gro -p ${outfolder}/topol.top -i ${outfolder}/posre.itp
        ## Move newly created files into "complex" folder

# Added line to fix topol_.itp files, to enable POSRES correctly
for f in `ls ${outfolder}/topol_*itp`; do line=$(grep "#include" $f); var=$(echo $line | awk -F'"' '{print $2}'); newvar=$(basename $var) ; newline=$(echo '#include "'${newvar}'"');  sed -i.bak "s~$line~$newline~g" $f  ; done


mv ${outfolder}/*itp ${outfolder}/complex/
mv ${outfolder}/topol.top ${outfolder}/complex/topol_pep.top
mv ${outfolder}/confout_pep.gro ${outfolder}/complex/

# Positioning of explicit bilayer at oriented protein position
	## Determine centre of implicit bilayer of oriented protein in [nm]
read x_imp y_imp z_imp < <( grep CA ${outfolder}/${protein}.pdb | awk '{x+=$7}{y+=$8}{z+=$9} END{print x/(10*NR),y/(10*NR),z/(10*NR)}') 

	## Determine centre of explicit bilayer in [nm]
	## Note1: Used P8 atoms only, with bilayer already in plane XY (Last frame)
	## Note2: Bilayer should be already solvated in electrolye, and equilibrated 
read x_exp y_exp z_exp < <(grep P8 ${DIR}/bilayer/mdf/confout_bil.pdb | awk '{x+=$6}{y+=$7}{z+=$8} END{print x/(10*NR),y/(10*NR),z/(10*NR)}')

	# Define translation vector for explicit bilayer
x_trans=$(echo "$x_imp - $x_exp" | bc -l)
y_trans=$(echo "$y_imp - $y_exp" | bc -l)
z_trans=$(echo "$z_imp - $z_exp" | bc -l)

	# Translate explicit bilayer centre (and whole box) to implicit bilayer centre
gmx_mpi editconf -f ${DIR}/bilayer/mdf/confout_bil.gro -o ${outfolder}/complex/confout_bil.gro -translate $x_trans $y_trans $z_trans

# Prepare GROMACS files for embedding with g_membed
cd complex
        # Get correct number of total atoms
pepfile=${outfolder}/complex/confout_pep.gro
Np="$(awk 'NR==2' $pepfile)"
bilfile=${outfolder}/complex/confout_bil.gro
Nb="$(awk 'NR==2' $bilfile)"
N=$((Np+Nb))

        ## Create and edit new file with protein and bilayer coordinates
echo "Protein in POPC in Water and Ions">> ${outfolder}/complex/for_embedding_NotCentred.gro
echo $N >> ${outfolder}/complex/for_embedding_NotCentred.gro
        ## Add protein coordinates, except headlines and last line
awk 'NR>2{print}' $pepfile | head -n -1 >> ${outfolder}/complex/for_embedding_NotCentred.gro
        ## Add bilayer(and electrolyte) coordinates, except headlines
        ## Preserve dimensions of solvated-bilayer box (last line)
awk 'NR>2{print}' $bilfile >> ${outfolder}/complex/for_embedding_NotCentred.gro
        ## Centre all coordinates in periodic box
gmx_mpi editconf -f ${outfolder}/complex/for_embedding_NotCentred.gro -o ${outfolder}/complex/for_embedding.gro -c
        ## Rename ions in .gro file, to be compatible with Berger f.f., naming
sed -i 's/K /K+/g' ${outfolder}/complex/for_embedding.gro
sed -i 's/CL /CL-/g' ${outfolder}/complex/for_embedding.gro
        ## Create index file with group "POPC_SOL_K+_CL-" and "SOL_K+_CL-" for embedding
printf "13|17|21|22\n17|21|22\nq\n" | gmx_mpi make_ndx -f ${outfolder}/complex/for_embedding.gro -o ${outfolder}/complex/for_embedding.ndx

