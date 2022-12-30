max_jobs=4; cur_jobs=0

for dirname in `ls -d data/cWza*_conformation*_[0-9][0-9][0-9][0-9]`; do
    echo $dirname;
    ((cur_jobs >= max_jobs)) && wait -n;
    time python md_contacts_analysis.py -v -p $dirname/md_100ns.tpr -x $dirname/md_100ns.xtc -o md_100ns.contacts.pickle & ((++cur_jobs));
done
