import numpy as np,healpy as hp,json
from scipy.interpolate import CubicSpline
from scipy.optimize import minimize
rng=np.random.default_rng(1)
T0=2.7255e6;beta=370/299792.458;gam=1/np.sqrt(1-beta**2)
psi=np.arange(0,181.0)
z=np.load(__import__('os').path.dirname(__file__)+'/../distributed/z4_dist_b_own.npz')['o'][:,1];T1=CubicSpline(psi,1/z)
PR=float(np.ptp(np.load(__import__('os').path.dirname(__file__)+'/../../../three_vel/point/cmb/obs_ring.npy')[:-1]))
def vec(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def lb(v):v=v/np.linalg.norm(v);return np.degrees(np.arctan2(v[1],v[0]))%360,np.degrees(np.arcsin(v[2]))
def ang(a,b):return np.degrees(np.arccos(np.clip(a@b/np.linalg.norm(a)/np.linalg.norm(b),-1,1)))
nside=32;V=np.array(hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))).T
aC=vec(0,15)  # anti-CoE in galactic; profile angle measured from anti-CoE
def sky(axis,v):  # observed T (relative) in galactic pixel frame; axis = anti-CoE unit vector, v = motion unit or None
    n=V
    if v is not None:
        mu=n@v;n=(n+((gam-1)*mu-gam*beta)[:,None]*v)/(gam*(1-beta*mu))[:,None];n/=np.linalg.norm(n,axis=1)[:,None]
        f=1/(gam*(1-beta*mu))
    else: f=1
    T=T1(np.degrees(np.arccos(np.clip(n@axis,-1,1))))*f
    return T/T.mean()-1
def dip(m):return 3*(V*m[:,None]).mean(0)       # dipole vector (dT/T)
def quad(m):
    a=hp.map2alm(m,lmax=2,iter=3);k=np.zeros_like(a)
    for mm in range(3):k[hp.Alm.getidx(2,2,mm)]=a[hp.Alm.getidx(2,2,mm)]
    q=hp.alm2map(k,nside);Q=np.einsum('p,pi,pj->ij',q,V,V)/len(q);e,E=np.linalg.eigh(Q)
    i=0 if abs(e[0]-e[1])>abs(e[2]-e[1]) else 2;return E[:,i],e
PH=vec(264.02,48.25);PQ=vec(257.73,57.4)
if __name__=="__main__":
    # roll check: axisymmetric by construction?
    m0=sky(aC,None);d0=dip(m0);q0,e0=quad(m0)
    out=dict(no_motion=dict(dipole_amp_uK=float(np.linalg.norm(d0)*T0),dipole_hot_lb=lb(d0),dipole_amp_over_Planck=float(np.linalg.norm(d0)*T0/3362),quad_axis_lb=lb(q0),quad_eigs=e0.tolist()))
    def fitv(axis,target):
        best=None
        for l0 in range(0,360,45):
            for b0 in (-60,0,60):
                r=minimize(lambda x:ang(dip(sky(axis,vec(*x))),target),[l0,b0],method='Nelder-Mead',options=dict(xatol=.1,fatol=1e-4))
                if best is None or r.fun<best.fun:best=r
        return best
    r=fitv(aC,PH);vf=vec(*r.x);m=sky(aC,vf);d=dip(m);q,e=quad(m)
    qo=ang(q,PQ);qo=min(qo,180-qo)
    out['fit']=dict(motion_lb=lb(vf),dipole_hot_lb=lb(d),dipole_offset_deg=float(r.fun),dipole_amp_uK=float(np.linalg.norm(d)*T0),dipole_amp_over_Planck=float(np.linalg.norm(d)*T0/3362),quad_axis_lb=lb(q),quad_offset_deg=float(qo),quad_axis_to_CoEaxis_deg=float(min(ang(q,aC),180-ang(q,aC))),range_over_Planck_ring=float(np.ptp(m)/PR))
    # also: motion-only dipole for reference
    out['motion_only_dipole_uK']=float(beta*T0)
    print(json.dumps(out,indent=1))
    # Monte Carlo: random isotropic CoE axis, motion fitted same way (quick fit: 1 start toward target)
    
    N=int(__import__('os').environ.get('N','500'));qs=[];ds=[]
    for i in range(N):
        a=rng.normal(size=3);a/=np.linalg.norm(a)
        rr=minimize(lambda x:ang(dip(sky(a,vec(*x))),PH),list(lb(PH)),method='Nelder-Mead',options=dict(xatol=.2,fatol=1e-3))
        mm=sky(a,vec(*rr.x));qq,_=quad(mm);o=ang(qq,PQ);qs.append(min(o,180-o));ds.append(rr.fun)
    qs=np.array(qs);p=(np.sum(qs<=qo)+1)/(N+1)
    from scipy.stats import norm
    out['MC']=dict(N=N,p=float(p),sigma_two_sided=float(norm.isf(p/2)),analytic_isotropic_p=float(1-np.cos(np.radians(qo))),median_dip_offset=float(np.median(ds)),frac_dip_offset_le_fit=float(np.mean(np.array(ds)<=r.fun+1e-6)))
    np.savez('mc.npz',q=qs,d=np.array(ds));json.dump(out,open('rotfit_stats.json','w'),indent=1,default=float);print(out['MC'])
    