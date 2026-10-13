import numpy as np,time
from run_asym import *
rng=np.random.default_rng(1); src=load_sources(); k='CCCCCCVVVVV'
p0=warm_p0(k,src[0],rng); m=Model(k); print('p0 left',m.leftover(p0),flush=True)
t=time.time(); m.pattern(p0); print('pattern s',time.time()-t,flush=True)
r=job((k,p0,'t')); print({a:r[a] for a in ['left','viol','asym_viol','mirror_matches','on_plane','closest_pair','min_plane_deg','odd_frac','secs']},flush=True)
