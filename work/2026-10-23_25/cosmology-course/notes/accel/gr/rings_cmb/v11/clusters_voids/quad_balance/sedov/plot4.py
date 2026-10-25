import os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
H=os.getcwd()
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
src=open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read()
exec(src.split('rng=np.random')[0]);os.chdir(H)
import sys;RUN=sys.argv[1];z0=np.load('runs/'+RUN+'.npz');rat=float(z0['ratio'])
rng=np.random.default_rng(12345);bb=None
for k in range(150):
    r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
    if bb is None or r.fun<bb.fun:bb=r
th=bb.x
for j in range(3): th=minimize(chi,th,args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-9,fatol=1e-9)).x
c2=chi(th,rat);m=model(th,rat);R2=1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)
# reference slc 0.041/0.035: direct refit
rr=1.063145729129294;rng=np.random.default_rng(12345);best=None
for k in range(40):
    r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rr,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
    if best is None or r.fun<best.fun:best=r
thr=best.x;c2r=chi(thr,rr)
Ls=np.mod(np.arange(0,361,1.0),360.0);Ls[-1]=360.0;Lsave=L.copy()
def curve(t,ra):
    global L;L=Ls;v=model(t,ra);L=Lsave;return v
def chi_at(t,ra,Lb):
    global L;L=np.mod(Lb,360.0);v=model(t,ra);L=Lsave;return np.sum(((v-Y)/S)**2)
print("bin L range",L.min(),L.max(),"chk combo",chi_at(th,rat,L),c2,"chk ref",chi_at(thr,rr,L),c2r)
from scipy.interpolate import interp1d
print("chi from plotted curve: combo",np.sum(((interp1d(Ls,curve(th,rat))(np.mod(L,360))-Y)/S)**2),"ref",np.sum(((interp1d(Ls,curve(thr,rr))(np.mod(L,360))-Y)/S)**2))
# CMB
z=np.load('runs/'+RUN+'.npz');o=z['o'];psi=o[:,0];lz=np.log(o[:,1])
c=np.cos(np.radians(psi));w=np.sin(np.radians(psi));dtt=-(lz-np.trapezoid(lz*w,psi)/np.trapezoid(w,psi))
from numpy.polynomial import legendre as Lg
from scipy.interpolate import CubicSpline;import healpy as hp,sys
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
mu=np.linspace(-1,1,4001);a=[0.0]
cs=CubicSpline(psi,lz);yy=-cs(np.degrees(np.arccos(mu)));yy-=0.5*np.trapezoid(yy,mu)
a+=[(2*l+1)/2*np.trapezoid(yy*Lg.legval(mu,[0]*l+[1]),mu) for l in range(1,7)]
T0=2.7255e6;ns=128;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T;pp=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));ll=np.arange(7)
def Dl(fn):
    o=np.load(fn)['o'];o=o[np.argsort(o[:,0])];c_=CubicSpline(o[:,0],np.log(o[:,1]));sm=np.exp(-c_(pp));sm=np.asarray(hp.remove_dipole(sm/sm.mean()-1))
    return ll*(ll+1)*hp.anafast(sm,lmax=6)*T0**2/(2*np.pi)
DC=Dl('runs/'+RUN+'.npz');DR_=Dl('../mass_family_scan/runs_sl2/slc_0.041037_0.035.npz')
PL=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt')[:5]
nodip=dtt-a[1]*c
fig,ax=plt.subplots(1,4,figsize=(25,5.5))
A=ax[0];A.errorbar(L,Y,yerr=S,fmt='o',ms=3,color='k',alpha=.6,label='Migkas 71 bins')
A.plot(Ls,curve(th,rat),'C3',lw=2,label=f'Sedov {RUN}: χ²={c2:.2f}, R²={R2:.3f}, ℓ0={th[1]:.1f}° (direct)')
A.plot(Ls,curve(thr,rr),'C0',ls='--',lw=1.6,zorder=5,label=f'ref slc 0.041/0.035: χ²={c2r:.2f}, ℓ0={thr[1]:.1f}° (direct)')
A.set_xlim(0,360);A.set_xlabel('galactic longitude l [deg]');A.set_ylabel('A/A_all');A.legend(fontsize=7);A.set_title('Migkas A/A_all')
A=ax[1];A.plot(psi,dtt*1e3,'C1o-',ms=3,label='ΔT/T incl. dipole');A.plot(psi,nodip*1e3,'C2s-',ms=3,label='dipole removed')
A.invert_xaxis();A.set_xlabel('ψ [deg] (180 = toward CoE, left)');A.set_ylabel('ΔT/T [10⁻³]')
A.set_title('CMB ΔT/T, 37 ray-traced 1+z');A.text(.03,.04,'\n'.join(f'a{l} = {a[l]:.3e}' for l in range(1,5)),transform=A.transAxes,fontsize=8);A.legend(fontsize=8)
x=np.linspace(1e-3,2.6,3000);Rb=1.5*(1+14.30048287260589/42)
exec(open('an.py').read().split('def shape')[1].split('rows=[]')[0].join(['def shape','']));sc=shape(z0,x);sr=x**-0.041037*np.exp(-0.035*(x-1)**2)*(x<=Rb)
A=ax[2];A.plot(x,sc,'C3',label='Sedov '+RUN+' ');A.plot(x,sr,'C0',ls='--',lw=1.6,zorder=5,label='ref slc 0.041/0.035')
A.axvline(1,color='k',ls=':');A.text(1.02,0.3,'observer',fontsize=8);A.plot(0,0,'k*',ms=12);A.text(0.03,0.05,'CoE',fontsize=8)
A.set_xlabel('comoving distance from CoE x (observer=1)');A.set_ylabel('initial ρ shape (ρ(1)≈1)');A.set_ylim(0,max(1.3,1.05*sc.max()));A.legend(fontsize=7,loc='upper left');A.set_title('Density profile')
I=inset_axes(A,width='40%',height='35%',loc='lower right',borderpad=2);m_=(x>0.85)&(x<1.15)
I.plot(x[m_],sc[m_],'C3');I.plot(x[m_],sr[m_],'C0',alpha=.35);I.axvline(1,color='k',ls=':');I.tick_params(labelsize=6)
A=ax[3];k=np.arange(2,7);A.bar(k-0.2,DC[2:],0.4,color='C3',label='Sedov '+RUN+'')
A.bar(k+0.2,DR_[2:],0.4,color='C0',alpha=.4,label='ref slc 0.041/0.035')
A.errorbar(k,PL[:,1],yerr=[PL[:,2],PL[:,3]],fmt='k_',ms=20,mew=2,capsize=4,label='Planck 2018 TT (R3.01 low-ℓ)')
A.set_yscale('log');A.set_xlabel('ℓ');A.set_ylabel('D_ℓ [μK²] (dipole removed)');A.legend(fontsize=7);A.set_title('Multipoles vs Planck')
plt.tight_layout();plt.savefig('best_4panel_'+RUN+'.png',dpi=110)
print(dict(chi2=c2,R2=R2,l0=th[1],ref_chi2=c2r,ref_l0=thr[1],a=a,dtt_range=[dtt.min(),dtt.max()],nodip_range=[nodip.min(),nodip.max()],Rb=Rb,DC=DC[2:].tolist(),DR=DR_[2:].tolist(),PL=PL[:,1].tolist()))
