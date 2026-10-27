# Diagnostic: power of compensation error (Resid - Tobs) and of body contributions split by distance from observer, l=10-40 (web_RA0.01_M1e22).
import numpy as np,healpy as hp,json
P5='/workspace/coepdf/joint_solve/phase5/';z=np.load(P5+'p5_web_RA0.01_M1e22.npz');y,Mb,B,Resid,Tobs,C=z['y'],z['M'],z['B'],z['Resid'],z['Tobs'],z['C']
NS=64;T0=2.7255e6;R=21.7439;EPS=0.025*R;MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
def kern(v):
    d=v[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@v/np.linalg.norm(v)**3;return (k-k.mean())*MSUN/KG*RV**2
def Dl(m,L=60):l=np.arange(L+1);return l*(l+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(m),lmax=L)*T0**2
r=np.linalg.norm(y,axis=1);out=dict(err=Dl(Resid-Tobs)[[15,20,24,30,35]].tolist(),resid=Dl(Resid)[[15,20,24,30,35]].tolist(),obs=Dl(Tobs)[[15,20,24,30,35]].tolist())
for lo,hi in [(0,.25),(.25,.5),(.5,.75),(.75,.875)]:
    s=np.where((r>=lo*R)&(r<hi*R))[0];m=sum(kern(y[i])*Mb[i] for i in s)
    out[f'bodies_r{lo}-{hi}R_n{len(s)}']=Dl(m)[[15,20,24,30,35]].tolist()
    out[f'resid_minus_r{lo}-{hi}R']=Dl(Resid-m)[[15,20,24,30,35]].tolist()
json.dump(out,open(P5+'web/ablation/diag_bins.json','w'),indent=1)
for k,v in out.items():print(k,np.round(v,1))
