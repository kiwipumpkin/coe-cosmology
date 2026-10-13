import numpy as np,sys
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast
for tag in sys.argv[1:]:
    f=np.load(f"final_{tag}.npz"); X0=f["Xp"]; r0=f["rp"]; M0=f["Mp"]
    dr=np.sqrt((np.hypot(X0[:,0],X0[:,1])-RMIN)**2+X0[:,2]**2); g=dr-r0
    print(tag,"SR+central alone: m13+ %.2f total %.1f"%(spec(R)["rms_13up"],spec(R)["total"]))
    for lo,hi in ((0,1),(1,3),(3,6),(6,10),(10,99)):
        m=(g>=lo)&(g<hi); W=(M0[m][:,None]*K(X0[m])).sum(0); sp=spec(W-W.mean())
        print("  bodies %4.0f-%2.0f Gpc from ring: n=%5d  m1-12 %.1f  m13+ %.2f"%(lo,hi,m.sum(),np.hypot(sp["rms_m1_5"],sp["rms_m6_12"]),sp["rms_13up"]))
