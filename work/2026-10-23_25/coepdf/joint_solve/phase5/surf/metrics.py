# Phase5 surface test metrics: 1+z range, map rms, dipole-removed rms, l=2 amplitude (sqrt(5 C2/4pi)... reported as rms of l=2 component), continuity (max |dx| between 5-deg rays).
import sys,json,numpy as np,healpy as hp
O='/workspace/coepdf/joint_solve/phase5/';R={}
for T in sys.argv[1:]:
    try:b=json.load(open(O+f'bg_{T}.json'));sm=np.load(O+f'bg_{T}_sm_ns64.npy');o=np.load(O+f'bg_{T}_trace.npz')['o']
    except Exception as e:R[T]=str(e);continue
    o=o[np.argsort(o[:,0])];a=hp.map2alm(sm,lmax=3,iter=3)
    m1=hp.alm2map(hp.almxfl(a,[0,1,0,0]),64);m2=hp.alm2map(hp.almxfl(a,[0,0,1,0]),64)
    dx=np.abs(np.diff(o[:,4]));R[T]=dict(surf=b.get('surf'),onepz=b['onepz_range'],R_ls=b['R_ls'],age=b['age_Gyr'],rms=float(sm.std()),rms_nodip=float((sm-m1).std()),l2_rms=float(m2.std()),
      x_hit=[float(o[:,4].min()),float(o[:,4].max())],max_dx=float(dx.max()),max_dx_at_psi=float(o[np.argmax(dx),0]),continuous=bool(dx.max()<0.1))
print(json.dumps(R,indent=1))
