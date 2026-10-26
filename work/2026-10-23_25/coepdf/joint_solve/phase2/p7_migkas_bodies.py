# Phase2: Migkas check with the phase-2 compensating bodies (same method as ../s5): dA/A ~ 2 dPhi/c^2 along each Migkas cone
# (b=0, |dl|<=45, 0.5-14.3 Gpc sightline average), applied to the locked eta-edge curve (d=14.3, D=34.96). Softening 0.1 Gpc as s5.
import sys,json,numpy as np
O='/workspace/coepdf/joint_solve/phase2/';G=4.30e-12;c2=(2.998e5)**2
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,yA=d['l_deg'],d['A_over_Aall'];s=0.5*(d['A_err_low']+d['A_err_high'])
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1);m0=np.interp(l,c[:,0],c[:,1],period=360)
rr=np.linspace(0.5,14.3,40);res={}
for f in sys.argv[1:]:
    Z=np.load(O+f);y,M=Z['y'],Z['M'];dPhi=[]
    for L0 in l:
        v=[]
        for dl in np.linspace(-45,45,19):
            t=np.radians(L0+dl);P=rr[:,None]*np.array([np.cos(t),np.sin(t),0.])
            v.append(np.mean(-G*np.sum(M[None]/np.sqrt(((P[:,None,:]-y[None])**2).sum(-1)+0.01),1)/c2))
        dPhi.append(np.mean(v))
    dPhi=np.array(dPhi);m1=m0*(1+2*(dPhi-dPhi.mean()))
    res[f]=dict(dPhi_over_c2_range=[float(dPhi.min()),float(dPhi.max())],chi2_before=float(np.sum((yA-m0)**2/s**2)),chi2_after=float(np.sum((yA-m1)**2/s**2)),
      r2_before=float(np.corrcoef(yA,m0)[0,1]**2),r2_after=float(np.corrcoef(yA,m1)[0,1]**2),max_dA=float(np.abs(m1-m0).max()))
    print(f,res[f],flush=True)
json.dump(res,open(O+'p7_migkas_bodies.json','w'),indent=1)
