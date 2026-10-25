import json,sys,numpy as np,healpy as hp
from core import *;from obs import *
NAMES=['C1','C2','C3','C4','log10frac','DRV','A_g'];L2=slice(2,LM+1);sig=np.sqrt(ClP[L2])
BASE={'red':dict(C=[-0.041,-0.035,0.1353,0],frac=20000,DRV=2.937,run='runs/red_base.npz'),
      'refpoly':dict(C=[-0.041,-0.035,0,0],frac=30,DRV=2.937,run='runs/refpoly.npz')}
STEP={'C1':('C',0,0.005),'C2':('C',1,0.01),'C3':('C',2,0.01),'C4':('C',3,0.01),'log10frac':('lf',0,0.1),'DRV':('DR',0,0.1)}
def model(run,C,DRV):
    o=np.load(run)['o'];sm,gm,_=terms(o,shape=shape_poly(C,DRV),DRV=DRV);return al0(sm),al0(gm)
Ga,Cl=obs_grid('/workspace/inv/smica_alm_l30.fits')
def bestaxis(am,ag,Ag=None):
    ao=Ga[L2];w=1/ClP[L2]
    if Ag is None:Ag=((ao-am[L2,None])*ag[L2,None]*w[:,None]).sum(0)/(ag[L2]**2*w).sum()
    c=((am[L2,None]+ag[L2,None]*Ag-ao)**2*w[:,None]).sum(0);i=int(np.argmin(c));return i,float(c[i]),float(np.broadcast_to(Ag,c.shape)[i])
if __name__=='__main__':
    R={}
    for U,b in BASE.items():
        am0,ag0=model(b['run'],b['C'],b['DRV']);Ag0=1.0;J=[]
        for n in NAMES[:-1]:
            k,i,h=STEP[n];C=list(b['C']);DRV=b['DRV']
            if k=='C':C[i]+=h
            if k=='DR':DRV+=h
            am,ag=model(f'runs/{U}_{"C%d"%(i+1) if k=="C" else k}.npz',C,DRV);J.append(((am+Ag0*ag)-(am0+Ag0*ag0))[L2]/h)
        J.append(ag0[L2]);J=np.array(J).T;Jw=J/sig[:,None]
        u,s,vt=np.linalg.svd(Jw,full_matrices=False);F=Jw.T@Jw
        cn=Jw/np.linalg.norm(Jw,axis=0);u2,s2,vt2=np.linalg.svd(cn,full_matrices=False)
        corr=cn.T@cn
        ia,c0,_=bestaxis(am0,ag0,1.0);ao=Ga[L2,ia];r=(ao-(am0+ag0)[L2])/sig
        steps={}
        for lam in [0,1e-3,1e-2,1e-1,1]:
            lamv=lam*s[0]**2
            d=np.linalg.solve(F+lamv*np.eye(len(NAMES)),Jw.T@r);steps[str(lam)]=d.tolist()
        R[U]=dict(names=NAMES,J_muK=J.tolist(),singular_values_weighted=s.tolist(),right_vectors_rows=vt.tolist(),
          singular_values_colnormalized=s2.tolist(),right_vectors_colnormalized=vt2.tolist(),column_norms=np.linalg.norm(Jw,axis=0).tolist(),
          column_cosine=corr.tolist(),rank_tol1em6=int((s>1e-6*s[0]).sum()),rank_tol1em3=int((s>1e-3*s[0]).sum()),
          cond=float(s[0]/s[-1]),start_axis=[float(x) for x in hp.pix2ang(NA,ia,lonlat=True)],start_chi2_Ag1=c0,LM_steps=steps,
          linear_pred_chi2={k:float(((r-Jw@np.array(v))**2).sum()) for k,v in steps.items()})
        json.dump(R,open('jac_results.json','w'),indent=1)
        print(U,'s',np.array2string(s,precision=3),'snorm',np.array2string(s2,precision=3),'start chi2',c0,R[U]['linear_pred_chi2'],flush=True)
        print(np.round(vt[-1],3),np.round(vt2[-1],3))
