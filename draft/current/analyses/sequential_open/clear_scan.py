"""Fine-scale texture vs a cleared zone around the ring: delete every body whose surface is within D of the ring; what m>=13 remains?"""
import numpy as np,sys,json
from common import *
from fit_on_web import spec
K=seq_core.unit_pattern_fast; out={}
for tag in ["b_random_s1_tidalcore","b_random_s2_tidalcore","b_random_s3_tidalcore","b_foam_rim_s1_tidalcore","b_foam_rim_s2_tidalcore","b_foam_rim_s3_tidalcore"]:
    f=np.load(f"final_{tag}.npz"); X0=f["Xp"]; r0=f["rp"]; M0=f["Mp"]; C=M0[:,None]*K(X0)
    g=np.sqrt((np.hypot(X0[:,0],X0[:,1])-RMIN)**2+X0[:,2]**2)-r0; out[tag]=[]
    for D in (0,1,2,3,4,5,6,7,8,10,12):
        keep=g>=D; res=R-C[keep].sum(0); sp=spec(res)
        out[tag].append(dict(D=D,deleted=int((~keep).sum()),frac=float((~keep).mean()),left=sp["total"],m1_12=float(np.hypot(sp["rms_m1_5"],sp["rms_m6_12"])),m13=sp["rms_13up"]))
        print(tag,D,out[tag][-1]["deleted"],"m13+ %.2f m1-12 %.0f"%(sp["rms_13up"],out[tag][-1]["m1_12"]))
json.dump(out,open("minedit2/clear_scan.json","w"),indent=1)
