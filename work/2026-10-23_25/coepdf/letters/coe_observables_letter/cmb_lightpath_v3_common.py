# Helper (no PNG): builds the v3A/B/C dT/T maps from the compute products and the shared statistics.
# Orientation: anti-CoE at (l,b)=(0,+15) (blended fig8 rows). All three fields are axisymmetric about the CoE-observer axis,
#   so the 80.8 deg roll and the gradient construction do not change them (they only re-label psi around the axis).
import numpy as np,healpy as hp,json
NS=256;T0uK=2.7255e6;LMAX=3*NS-1
anti=np.array(hp.ang2vec(0.,15.,lonlat=True));nv=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T
psi_pix=np.degrees(np.arccos(np.clip(nv@anti,-1,1)))
def profile(rule):
    """returns psi_deg, ratio(psi) = T_obs/T_rec (A,B) or linear dT/T (C), and kind"""
    if rule=='A':
        z=np.load('cmb_lightpath_v3A.npz');onez=-z['pt_src'];return z['psi_deg'],1/onez,'ratio'
    z=np.load('cmb_lightpath_v3BC.npz')
    if rule=='B':
        r=z['resB'];return r[:,0],-r[:,5]+r[:,6]+r[:,7],'log'   # ln(T/T_rec); grav+ISW ~1e7 so exp() overflows: map = ln T - <ln T> (first order)
    r=z['resC'];return r[:,0],r[:,5]+r[:,6]+r[:,7],'linear'
def total_map(rule):
    p,y,kind=profile(rule);ok=np.isfinite(y);m=np.interp(psi_pix,p[ok],y[ok])
    if kind=='ratio':return m/m.mean()-1,kind,float(m.mean())
    return m-m.mean(),kind,float(m.mean())   # 'linear' (C) and 'log' (B)
def stats(rule):
    m,kind,mean=total_map(rule);m1=hp.remove_dipole(m)
    cl=hp.anafast(m,lmax=LMAX,iter=3);cl1=hp.anafast(m1,lmax=LMAX,iter=3)
    return m,m1,cl,cl1,kind,mean
def planck():
    p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');return p[:,0],p[:,1]
def pct(cl1):
    l=np.arange(2,31);c=cl1[2:31]*(2*l+1);L,Dl=planck();s=(L>=2)&(L<=30);P=(2*L[s]+1)*2*np.pi*Dl[s]/(L[s]*(L[s]+1))
    return l,100*c/c.sum(),100*P/P.sum()
def axis(m,l):
    alm=hp.map2alm(m,lmax=l,iter=3);a=np.zeros(2*l+1,complex)
    for mm in range(l+1):
        v=alm[hp.Alm.getidx(l,l,mm)];a[l+mm]=v
        if mm:a[l-mm]=(-1)**mm*np.conj(v)
    ms=np.arange(-l,l+1);Lz=np.diag(ms).astype(complex);Lp=np.zeros((2*l+1,)*2,complex)
    for i,mm in enumerate(ms[:-1]):Lp[i+1,i]=np.sqrt(l*(l+1)-mm*(mm+1))
    Lx=(Lp+Lp.conj().T)/2;Ly=(Lp-Lp.conj().T)/2j;Ls=[Lx,Ly,Lz]
    A=np.array([[np.real(np.vdot(Li@a,Lj@a)) for Lj in Ls] for Li in Ls]);w,v=np.linalg.eigh(A);n=v[:,-1]
    lo,la=hp.vec2ang(n,lonlat=True);aoe=np.array(hp.ang2vec(260,60,lonlat=True))
    return dict(l=float(lo[0]),b=float(la[0]),angle_to_AoE=float(np.degrees(np.arccos(abs(n@aoe)))),angle_to_CoE_axis=float(np.degrees(np.arccos(abs(n@anti)))),
                eig=w.tolist(),degenerate=bool(abs(w[-1]-w[-2])<1e-6*abs(w[-1])))
