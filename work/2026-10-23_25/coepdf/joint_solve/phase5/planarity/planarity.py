# Planarity / alignment test of low-ell CMB (Rakic&Schwarz 2007; de Oliveira-Costa et al 2004) for Planck target, v3 seed-1 model sky, CoE-only bg.
# Skies from ../web/inverse/v3/inv3_s1.npz: Tobs (SMICA l=2..130, Nside64, full sky no mask), Resid (bg+web final), B (bg only), R0 (bg+random web start).
import numpy as np,healpy as hp,json
V='/workspace/coepdf/joint_solve/phase5/web/inverse/v3/';T0=2.7255e6;LM=5;AOE=(260.,60.);NSIM=10000
d=np.load(V+'inv3_s1.npz');sk={'Planck_target':d['Tobs'],'v3_s1_model':d['Resid'],'CoE_bg_only':d['B'],'v3_s1_start':d['R0']}
def full(alm,l):  # complex a_lm, m=-l..l
    p=np.array([alm[hp.Alm.getidx(LM,l,m)] for m in range(l+1)]);return np.r_[[(-1)**m*np.conj(p[m]) for m in range(l,0,-1)],p]
def Lmat(l):
    m=np.arange(-l,l+1);Lz=np.diag(m).astype(complex);Lp=np.zeros((2*l+1,)*2,complex)
    for i,mm in enumerate(m[:-1]):Lp[i+1,i]=np.sqrt(l*(l+1)-mm*(mm+1))
    Lm=Lp.conj().T;return [(Lp+Lm)/2,(Lp-Lm)/(2j),Lz]
LL={l:Lmat(l) for l in range(2,LM+1)}
def bestaxis(v,l):  # maximise n.M.n, M_ij = Re v^H L_i L_j v (exact eigenproblem)
    L=LL[l];M=np.array([[np.real(v.conj()@L[i]@L[j]@v) for j in range(3)] for i in range(3)]);M=(M+M.T)/2
    w,e=np.linalg.eigh(M);n=e[:,-1];return n if n[2]>=0 else -n,w[-1]/np.real(v.conj()@v)
def rot_to(alm,lon,lat):  # rotate so that (lon,lat) -> north pole; verified below
    a=alm.copy();hp.rotate_alm(a,-np.radians(lon),-np.radians(90-lat),0.,lmax=LM) if False else None
    r=hp.Rotator(rot=[lon,lat-90.],deg=True,eulertype='ZYX');return r.rotate_alm(alm,lmax=LM)
# verify rotation convention with a spot
spot=np.zeros(hp.nside2npix(64));spot[hp.ang2pix(64,*AOE,lonlat=True)]=1;sa=hp.map2alm(hp.smoothing(spot,fwhm=0.3),lmax=LM)
pk=hp.pix2ang(64,np.argmax(hp.alm2map(rot_to(sa,*AOE),64,lmax=LM)),lonlat=True);assert pk[1]>80,('rotation check failed',pk)
def ang(u,v):return float(np.degrees(np.arccos(min(1,abs(u@v)))))
naoe=np.array(hp.ang2vec(*AOE,lonlat=True))
def stats(alm):
    ar=rot_to(alm,*AOE);o={};ax={}
    for l in range(2,LM+1):
        v=full(ar,l);p=np.abs(v)**2;S=p.sum();n,Lz2=bestaxis(full(alm,l),l);ax[l]=n;lo,la=hp.vec2ang(n,lonlat=True)
        o[l]=dict(t_l=float((p[0]+p[-1])/S),zonal=float(p[l]/S),C_l_muK2=float(S/(2*l+1)*T0**2),best_axis_lb=[float(lo[0]),float(la[0])],Lz2_norm=float(Lz2),ang_to_AoE=ang(n,naoe))
    o['ang_l2_l3']=ang(ax[2],ax[3]);return o,ax
R={k:stats(hp.map2alm(m,lmax=LM,iter=3))[0] for k,m in sk.items()}
# isotropic Gaussian sims: t_l about fixed axis (frame-independent), l2-l3 best-axis angle, Lz2_norm
rng=np.random.default_rng(0);sim={l:[] for l in range(2,LM+1)};sang=[];slz={l:[] for l in range(2,LM+1)}
for i in range(NSIM):
    ax={}
    for l in range(2,LM+1):
        p=rng.normal(size=l+1)+1j*rng.normal(size=l+1);p[0]=p[0].real*np.sqrt(2)
        v=np.r_[[(-1)**m*np.conj(p[m]) for m in range(l,0,-1)],p];q=np.abs(v)**2;sim[l].append((q[0]+q[-1])/q.sum())
        if l<=3:n,lz=bestaxis(v,l);ax[l]=n;slz[l].append(lz)
    sang.append(ang(ax[2],ax[3]))
sang=np.array(sang)
for k in R:
    for l in range(2,LM+1):
        R[k][l]['p_t_l_ge']=float(np.mean(np.array(sim[l])>=R[k][l]['t_l']))
        if l<=3:R[k][l]['p_Lz2_ge']=float(np.mean(np.array(slz[l])>=R[k][l]['Lz2_norm']))
    R[k]['p_ang_l2_l3_le']=float(np.mean(sang<=R[k]['ang_l2_l3']))
R['sim_note']=f'{NSIM} isotropic Gaussian skies; t_l and alignment are C_l-independent so C_l scaling irrelevant. Full sky, no mask.'
json.dump(R,open('planarity.json','w'),indent=1);print(json.dumps(R,indent=1))
