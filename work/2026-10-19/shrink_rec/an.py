import numpy as np,json,glob
LD=18159.848409348455/39.10676196793444;al=0.30000009809958916;q=(1/LD)**al;nlaw=al*q/(1+q)
def dip(o):
    z=o[:,1];p=np.radians(o[:,0]);w=np.sin(p);w[0]=w[-1]=np.sin(np.radians(2.5))/4;w/=w.sum()
    T=1/z;d=T/np.sum(w*T)-1;X=np.c_[np.ones_like(p),np.cos(p)];sw=np.sqrt(w);c=np.linalg.lstsq(X*sw[:,None],d*sw,rcond=None)[0];return c[1],1/np.sum(w*T)
R={}
for f in sorted(glob.glob('*_k*.npz')):
    Z=np.load(f);o=Z['o'];c,m=dip(o);te=o[:,2];XX,tt,Rr=Z['XX'],Z['te'],Z['Rr']
    i=(XX>0.85)&(XX<1.2);sl=np.polyfit(tt[i],Rr[i],1)[0]   # dR_rec/dt across visible shells
    R[f[:-4]]=dict(slope_n=float(nlaw+float(Z['k'])),xPlanck=float(abs(c)*2.7255e6/3362),hot_at_CoE=bool(c<0),mean_1pz=float(m),te_CoE=float(te[-1]),te_antiCoE=float(te[0]),te_spread_pct=float((te[-1]-te[0])/te.mean()*100),dRrec_dt=float(sl),tus=float(Z['tus']))
    print(f,{k:round(v,5) if isinstance(v,float) else v for k,v in R[f[:-4]].items()})
print('n_law at x=1',nlaw);json.dump(dict(runs=R,n_law=nlaw),open('shrink_results.json','w'),indent=1)
