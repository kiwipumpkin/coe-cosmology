import numpy as np,sys
def dip(f):
    o=np.load(f)['o'];ok=np.isfinite(o[:,1]);z=o[:,1];p=np.radians(o[:,0]);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(.5))/4;w/=w.sum()
    T=1/z;d=T/np.sum(w*T)-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0]
    return dict(n_ok=int(ok.sum()),mean_1pz=float(1/np.sum(w*T)),dip=float(c[1]),xPlanck=float(abs(c[1])*2.7255e6/3362),te_min=float(np.nanmin(o[:,2])),te_max=float(np.nanmax(o[:,2])))
if __name__=='__main__':
    for f in sys.argv[1:]:print(f,dip(f))
