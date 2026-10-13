cd /workspace/clump_scan/inside_visible
NP=6 python3 fit_in.py cap2e52 40 2e52 $(cat classes.txt) 2>&1 | grep -v beta
NP=6 python3 fit_in.py cap1e53 40 1e53 $(cat classes.txt | tr ' ' '\n' | grep C | tr '\n' ' ') 2>&1 | grep -v beta
MASS="C CC CCC CCCC CCCCC CCCCCC CCCCCCC CCCCCCCC CCCCCCCCC CCCCCCCCCC"
COMP=none NP=6 python3 fit_in.py nocomp 40 2e52 $MASS 2>&1 | grep -v beta
COMP=surround NP=6 python3 fit_in.py surrcomp 40 2e52 $MASS 2>&1 | grep -v beta
