cd /workspace/clump_scan/inside_visible
while kill -0 1366373 2>/dev/null; do sleep 10; done
MASS="C CC CCC CCCC CCCCC CCCCCC CCCCCCC CCCCCCCC CCCCCCCCC CCCCCCCCCC"
NP=4 python3 -u fit_in.py cap1e53 12 1e53 $MASS CCCCV CCCCVV
COMP=none NP=4 python3 -u fit_in.py nocomp 12 2e52 C CC CCC CCCC CCCCC CCCCCC CCCCCCC CCCCCCCC
COMP=surround NP=4 python3 -u fit_in.py surrcomp 12 2e52 C CC CCC CCCC CCCCC CCCCCC CCCCCCC CCCCCCCC
