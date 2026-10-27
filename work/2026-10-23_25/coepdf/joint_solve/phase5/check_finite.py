# Phase5 (NEW): exit 1 unless every float array in each given .npz/.npy exists and is finite. Used by cmb_mig.sh to fail loudly.
import sys,numpy as np,os
bad=0
for f in sys.argv[1:]:
    if not os.path.exists(f):print('MISSING',f);bad=1;continue
    z=np.load(f);it=z.items() if hasattr(z,'files') else [('arr',z)]
    for k,v in it:
        if v.dtype.kind in 'fc' and not np.all(np.isfinite(v)):print('NONFINITE',f,k,int((~np.isfinite(v)).sum()));bad=1
sys.exit(bad)
