# usage: poly_plots.py PASS.npz BESTCHI.npz   (both smooth-polynomial runs)
import os,sys,json,numpy as np,warnings;warnings.filterwarnings('ignore')
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.getcwd();FP,FB=[os.path.abspath(a) for a in sys.argv[1:3]];OUT=sys.argv[3];REF=os.path.abspath('../runs_sl2/slc_0.041037_0.035.npz')
os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
src=open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read()
exec(src.split('rng=np.random')[0]);os.chdir(H)
sys.argv=[sys.argv[0]];import anP,glob
def mfit(rat,n=40):
    rng=np.random.default_rng(12345);best=None
    for k in range(n):
        r=minimize(chi,LO+(HI-LO)*rng.random(6),args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-7,fatol=1e-7))
        if best is None or r.fun<best.fun:best=r
    th=best.x
    for j in range(3):th=minimize(chi,th,args=(rat,),method='Nelder-Mead',options=dict(maxiter=5000,xatol=1e-9,fatol=1e-9)).x
    m=model(th,rat);return th,chi(th,rat),1-np.sum((Y-m)**2)/np.sum((Y-Y.mean())**2)
Ls=np.arange(0,361,1.0);Lsave=L.copy()
def curve(t,ra):
    global L;L=np.mod(Ls,360);L[-1]=360.0;v=model(t,ra);L=Lsave;return v
from scipy.interpolate import interp1d
out={}
for tag,fn in [('pass',FP),('bestchi',FB),('ref',REF)]:
    A=anP.an(fn);d0=np.load(fn);A['osc']=[float(d0[k]) for k in ('A','LAM','LD','PH')] if 'A' in d0.files else None;th,c2,R2=mfit(A['ratio']);cc=float(np.sum(((interp1d(Ls,curve(th,A['ratio']))(np.mod(L,360))-Y)/S)**2))
    d=np.load(fn);A.update(chi2=float(c2),chi2_from_plotted=cc,R2=float(R2),l0=float(th[1]),theta=th.tolist(),frac=float(d['frac']));out[tag]=A
    print(tag,A['C'],'ratio %.5f chi2 %.3f (plot %.3f) R2 %.3f l0 %.2f a'%(A['ratio'],c2,cc,R2,th[1]),['%.3e'%x for x in A['a']],'D2-6',['%.3g'%x for x in A['Dl'][2:7]])
PL=np.loadtxt('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb/planckTT.txt')
out['planck_l2_30']=PL[:29].tolist();scan=[]
for f in sorted(glob.glob('runs_osc/*.npz')):
    r=anP.an(f);r['tag']=os.path.basename(f)[:-4];r.pop('Dl_full',None);scan.append(r)
out['scan']=scan;json.dump(out,open(OUT+'_results.json','w'),indent=1)
col={'pass':'C3','bestchi':'C2','ref':'C0'};lab={k:(os.path.basename(out[k]['run'])[:-4] if k!='ref' else 'ref slc 0.041/0.035')+('') for k in out if k in col}
fig,ax=plt.subplots(1,4,figsize=(25,5.5))
A=ax[0];A.errorbar(L,Y,yerr=S,fmt='o',ms=3,color='k',alpha=.6,label='Migkas 71 bins')
for k in col:A.plot(Ls,curve(np.array(out[k]['theta']),out[k]['ratio']),col[k],ls='-' if k=='pass' else '--',lw=1.8,label=lab[k]+f": χ²={out[k]['chi2']:.2f} (plot {out[k]['chi2_from_plotted']:.2f}), ℓ0={out[k]['l0']:.1f}°")
A.set_xlim(0,360);A.set_xlabel('galactic longitude l [deg]');A.set_ylabel('A/A_all');A.legend(fontsize=6);A.set_title('Migkas A/A_all (direct fits)')
A=ax[1]
for k in col:
    o=np.load(out[k]['run'])['o'];o=o[np.argsort(o[:,0])];psi=o[:,0];lz=np.log(o[:,1]);w=np.sin(np.radians(psi));d_=-(lz-np.trapezoid(lz*w,psi)/np.trapezoid(w,psi))
    A.plot(psi,d_*1e3,col[k],marker='o',ms=2,label=lab[k])
A.invert_xaxis();A.set_xlabel('ψ [deg] (180 = toward CoE, left)');A.set_ylabel('ΔT/T [10⁻³] (monopole removed)');A.set_title('CMB ΔT/T, ray-traced')
p=out['pass'];A.text(.03,.04,'pass: '+', '.join(f'a{l+1}={p["a"][l]:.2e}' for l in range(4)),transform=A.transAxes,fontsize=7);A.legend(fontsize=6)
A=ax[2];x=np.linspace(1e-3,2.6,3000);Rb=1.5*(1+14.30048287260589/42)
for k in col:
    import re;tg=os.path.basename(out[k]['run']);m_=re.match(r'o_([^_]+)_([^_]+)_([^_]+)_([a-z]+)',tg)
    if m_:
        a_,lm,ld,pt=float(m_[1]),float(m_[2]),float(m_[3]),m_[4];ph={'trough':-np.pi/2,'crest':np.pi/2,'up':0,'down':np.pi}[pt]-2*np.pi/lm;xx=np.minimum(x,Rb);s=xx**-0.041037*np.exp(-0.035*(xx-1)**2+a_*np.exp(-xx/ld)*np.sin(2*np.pi*xx/lm+ph))
    else:s= x**-0.041037*np.exp(-0.035*(x-1)**2)
    A.plot(x,s*(x<=Rb),col[k],label=lab[k])
A.axvline(1,color='k',ls=':');A.text(1.02,0.1,'observer',fontsize=8);A.plot(0,0,'k*',ms=12);A.text(0.03,0.05,'CoE',fontsize=8)
A.set_xlabel('comoving distance from CoE x (observer=1)');A.set_ylabel('initial ρ shape (ρ(1)=1)');A.legend(fontsize=6,loc='upper left');A.set_title('Density vs x (CoE left)')
A=ax[3];k_=np.arange(2,7)
for i,k in enumerate(col):A.bar(k_-0.27+0.27*i,np.array(out[k]['Dl'][2:7]),0.27,color=col[k],alpha=.8,label=lab[k])
A.errorbar(k_,PL[:5,1],yerr=[PL[:5,2],PL[:5,3]],fmt='k_',ms=20,mew=2,capsize=4,label='Planck 2018 TT')
A.set_yscale('log');A.set_xlabel('ℓ');A.set_ylabel('D_ℓ [μK²] (dipole removed)');A.legend(fontsize=6);A.set_title('Multipoles vs Planck')
plt.tight_layout();plt.savefig(OUT+'_plots.png',dpi=90)
plt.figure(figsize=(9,5.5));ll=np.arange(2,31)
for k in col:plt.semilogy(ll,np.maximum(np.array(out[k]['Dl'][2:31]),1e-6),col[k],marker='o',ms=3,label=lab[k])
plt.errorbar(PL[:29,0],PL[:29,1],yerr=[PL[:29,2],PL[:29,3]],fmt='ks',ms=3,capsize=2,label='Planck 2018 TT')
plt.xlabel('ℓ');plt.ylabel('D_ℓ [μK²]');plt.title('Angular power spectrum ℓ=2–30 (37 rays, 5° spacing; ℓ≳20 resolution-limited)');plt.legend(fontsize=7);plt.grid(alpha=.3)
plt.tight_layout();plt.savefig(OUT+'_spectrum.png',dpi=90)
