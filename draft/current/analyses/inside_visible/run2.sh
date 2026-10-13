cd /workspace/clump_scan/inside_visible
while kill -0 1336980 2>/dev/null; do sleep 10; done
MASS="C CC CCC CCCC CCCCC CCCCCC CCCCCCC CCCCCCCC CCCCCCCCC CCCCCCCCCC"
NP=6 python3 -u fit_in.py cap1e53 20 1e53 $MASS CCCCV CCCCVV CCCCCVV CCCCVVVV
COMP=none NP=6 python3 -u fit_in.py nocomp 20 2e52 $MASS
COMP=surround NP=6 python3 -u fit_in.py surrcomp 20 2e52 $MASS
