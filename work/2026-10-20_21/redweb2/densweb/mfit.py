import sys,json,csv,numpy as np,healpy as hp
from scipy.optimize import minimize_scalar
W_=np.load('web_v_s1.npz');DG=W_['DG'];Vp=W_['Vp'];Vh=W_['Vh']
DG=np.r_[DG,2.5,3.0];Vp=np.pad(Vp,((0,0),(0,2)),mode='edge');Vh=np.pad(Vh,((0,0),(0,2)),mode='edge')  # web held at 2 Gpc value beyond grid (only z>0.40 clusters)
T=np.load('ltb_t0.npz');xg,Rg,Rtg=T['x'],T['R'],T['Rt'];Rt1=np.interp(1,xg,Rtg);H1=Rt1/np.interp(1,xg,Rg)
DM=42000.;C=299792.458;SC=70/H1*DM
vec=lambda l,b:np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))*np.ones_like(np.asarray(l,float))])
a=vec(0.,15.)           # direction to CoE (run11: aG=fit.vec(0,15))
def cz_smooth(N):       # N (n,3) unit sightlines; returns (n,len(DG)) km/s
    O=-a;out=np.empty((len(N),len(DG)))
    for k,d in enumerate(DG/42.):
        P=O[None]+d*N;r=np.linalg.norm(P,axis=1);x=np.interp(r,Rg,xg);v=np.interp(x,xg,Rtg)[:,None]*P/r[:,None]
        out[:,k]=((v-Rt1*O[None])*N).sum(1)*SC
    return out
M=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
L=np.array([float(r['l_deg']) for r in M]);Y=np.array([float(r['A_over_Aall']) for r in M]);S=0.5*(np.array([float(r['A_err_low']) for r in M])+np.array([float(r['A_err_high']) for r in M]));Wt=1/S**2
Tc=list(csv.DictReader(open('/workspace/cosmology/migkas/migkas2020_tableC1.csv')));CZ=C*np.array([float(r['z']) for r in Tc]);CL=np.array([float(r['l']) for r in Tc]);CB=np.array([float(r['b']) for r in Tc])
BINS=[np.where(np.abs((CL-x+180)%360-180)<=45)[0] for x in L]
chi=lambda m:float(np.sum(Wt*(m-Y)**2));r2=lambda m:float(np.corrcoef(m,Y)[0,1]**2) if np.std(m)>0 else 0.
def place(CZt):  # CZt (ndir,nd) -> H(ndir,ncl) km/s/Mpc
    H=np.empty((len(CZt),len(CZ)))
    for j,c in enumerate(CZt):
        if np.any(np.diff(c)<=0) or c[0]>CZ.min() or c[-1]<CZ.max():raise RuntimeError('non-monotonic/unbracketed cz on dir %d'%j)
        H[j]=CZ/(np.exp(np.interp(np.log(CZ),np.log(c),np.log(DG)))*1e3)
    return H
NP=vec(np.arange(360.),0*np.arange(360.)).T;SP=cz_smooth(NP)
NH=np.array(hp.pix2vec(16,np.arange(hp.nside2npix(16)))).T;SH=cz_smooth(NH)
def plane(s):
    G=place(SP+s*Vp);return G/G.mean(0)
def mplane(G,D):
    out=np.empty(len(L))
    for b in range(len(L)):
        t=(L[b]-D)%360;i0=int(np.floor(t))%360;f=t-np.floor(t);out[b]=np.mean((1-f)*G[i0,BINS[b]]+f*G[(i0+1)%360,BINS[b]])
    return out
def mown(s,D):
    Gn=place(SH+s*Vh).mean(0)              # full-sky mean per cluster
    th=np.radians(90-CB);ph=np.radians(CL-D)
    cz=np.array([np.interp(0,[0],[0])]) # dummy
    smo=cz_smooth(vec(CL-D,CB).T)
    web=np.array([hp.get_interp_val(Vh[:,k],th,ph) for k in range(len(DG))]).T
    H=np.array([CZ[i]/(np.exp(np.interp(np.log(CZ[i]),np.log(smo[i]+s*web[i]),np.log(DG)))*1e3) for i in range(len(CZ))])/Gn
    return np.array([H[ix].mean() for ix in BINS])
def hump(m):
    o=np.argsort(L);ll,mm=L[o]%360,m[o];pk=[float(ll[i]) for i in range(len(mm)) if mm[i]>=mm[i-1] and mm[i]>=mm[(i+1)%len(mm)]]
    return pk,float(ll[np.argmax(mm)]),float(ll[np.argmin(mm)])
res={}
for tag,s in [('smooth',0.),('web1.5',.015),('web1.75',.0175),('web2',.02)]:
    G=plane(s);sc=[chi(mplane(G,D)) for D in range(360)];D0=float(np.argmin(sc))
    r=minimize_scalar(lambda D:chi(mplane(G,D)),bounds=(D0-1,D0+1),method='bounded');Df=float(r.x%360)
    mf=mplane(G,Df);m0=mplane(G,0.)
    so=[chi(mown(s,D)) for D in range(0,360,3)];Do=3*float(np.argmin(so));ro=minimize_scalar(lambda D:chi(mown(s,D)),bounds=(Do-3,Do+3),method='bounded');mo=mown(s,ro.x);mo0=mown(s,0.)
    d=dict(strength=s,
      plane_fit=dict(CoE_l=Df,l0_lowz_convention=(Df+180)%360,chi2=chi(mf),pearson_r2=r2(mf),amp_half_ptp=float(np.ptp(mf)/2),humps=hump(mf),A=mf.tolist()),
      plane_own_axis=dict(CoE_l=0.,chi2=chi(m0),pearson_r2=r2(m0),amp_half_ptp=float(np.ptp(m0)/2),humps=hump(m0),A=m0.tolist()),
      owndir_fit=dict(CoE_l=float(ro.x%360),chi2=chi(mo),pearson_r2=r2(mo),amp_half_ptp=float(np.ptp(mo)/2),humps=hump(mo),A=mo.tolist()),
      owndir_own_axis=dict(CoE_l=0.,chi2=chi(mo0),pearson_r2=r2(mo0),amp_half_ptp=float(np.ptp(mo0)/2),A=mo0.tolist()),
      max_web_vlos_kms_le1Gpc=float(np.abs(s*Vp[:,DG<=1]).max()),smooth_cz_aniso_kms_at_0p3Gpc=float(np.ptp(SP[:,np.argmin(abs(DG-.3))])))
    res[tag]=d;print(tag,{k:(v['chi2'],v['pearson_r2'],v['amp_half_ptp'],v['CoE_l']) for k,v in d.items() if isinstance(v,dict)},d['max_web_vlos_kms_le1Gpc'],flush=True)
res['flat_chi2']=chi(np.ones_like(Y));res['notes']='RED LTB R(t0,x) from famP poly path (ratio %.7f reproduced); positions/velocities at t0, lowz-style; web v_pec=2f g/(3H Om), masses x strength'%float(T['ratio'])
json.dump(res,open('/workspace/redweb2/densweb/migkas_densweb_s1.json','w'),indent=1)
