
n_cores=6
for dir in `ls -d data/cWza*_conformation[0,1]_[0-9][0-9][0-9][0-9]`; do
	FILE=${dir}/frames
	if [ -d "$FILE" ]; then
	   echo "$FILE exists."
	else
	   echo "INFO: Extracting file in parallel."
	   python extract_frames_parallel.py $dir $n_cores
        fi
	zip -r "${FILE}".zip $FILE
	rm -rf "${FILE}"/frame*.pdb
done
