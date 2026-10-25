import sys,glob,json,numpy as np;sys.path.insert(0,'../../dipole_amp')
def dip(f):
    o=np.load(f)['o'];ok=np.isfinite(o[:,1]);o=o[ok];z=o[:,1];p=np.radians(o[:,0]);w=np.sin(p);w[w==0]=np.sin(np.radians(2.5))/4;w/=w.sum()
    T=1/z;d=T/np.sum(w*T)-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0]
    return dict(n_ok=int(ok.sum()),mean_1pz=float(1/np.sum(w*T)),dip=float(c[1]),xPlanck=float(abs(c[1])*2.7255e6/3362),hot_at_CoE=bool(c[1]<0),w_antiCoE=float(o[0,4]),w_CoE=float(o[-1,4]),te=float(o[0,2]))
R={f[4:-4]:dip(f) for f in sorted(glob.glob('sim_*.npz'))};R['shrink_b_1100']=dip('../b_1100_k0.npz')
for k,v in R.items():print(k,{a:round(b,5) if isinstance(b,float) else b for a,b in v.items()})
json.dump(R,open('floor_results.json','w'),indent=1)
