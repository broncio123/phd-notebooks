for dirname in `ls -d data/cWza*_conformation*_[0-9][0-9][0-9][0-9]`; do 
    echo $dirname;
    time python md_hbonds_analysis.py -v -p $dirname/md_100ns.tpr -x $dirname/md_100ns.xtc -o md_100ns.hbonds.pickle
done
