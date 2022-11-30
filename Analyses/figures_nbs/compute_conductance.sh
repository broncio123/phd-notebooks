n_cpus=7
for dirname in `ls -d data/cWza*_conformation*_[0-9][0-9][0-9][0-9]`; do 
    echo $dirname;
    rm -rf $dirname/frames;
    unzip -j $dirname/frames.zip -d $dirname/frames;
    time ls $dirname/frames/*pdb | python run_hole_parallel.py $n_cpus $dirname;
    rm -rf $dirname/frames;
done
