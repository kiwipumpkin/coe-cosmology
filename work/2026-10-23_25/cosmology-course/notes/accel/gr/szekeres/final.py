import pickle,json,glob,sys,numpy as np,warnings; warnings.filterwarnings('ignore'); sys.path.insert(0,'..'); sys.path.insert(0,'/workspace/cosmology-course/notes/paper_repro/grfit')
from grcore import stats,norm
from core import TH1,IMIG,Y,MLO,MHI
import healpy as hp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
T0=2.7255; u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dv=u(264.021,48.253); beta=3362.08e-6/T0; Pp=dv-(dv@A)*A; Pp/=np.linalg.norm(Pp)
th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),Pp)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); Mk=hp.read_map('/workspace/planck/mask.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.empty(361); mf=np.empty(361)
for i,v in enumerate(n): pix=hp.query_disc(ns,v,np.radians(1.0),nest=True); cmb[i]=I[pix].mean(); mf[i]=1-Mk[pix].mean()
bad=(mf>0.5)|(np.abs(np.degrees(np.arcsin(n[:,2])))<10); full=cmb/T0+beta*(n@dv); full-=full[:-1].mean()
rms=lambda x,m=None: np.sqrt(np.mean((x[:-1] if m is None else x[:-1][~m[:-1]])**2))
M=pickle.load(open('mig.pkl','rb')); R={}
for f in sorted(glob.glob('cmb_*.pkl')):
    d=pickle.load(open(f,'rb')); x=norm(d['x']); s=stats(x); m=M[(d['i'],d['C'])]
    R[(d['i'],d['C'])]=dict(x=x,dip=s[0],quad=s[1],cold=x.min(),rms=rms(x),left=rms(x-full),leftu=rms(x-full,bad),chi2=m['chi2'],R2=m['R2'],af=m['af'],sig0=d['sig0'],ratio=d['ratio'])
for k in sorted(R): r=R[k]; print(k,'max|S\'/S| %.4f (%.0f%% of no-crossing limit)'%(abs(k[1])*r['sig0'],100*r['ratio']),'dip %.5f quad %.5f cold %.5f rms %.5f left %.5f/%.5f chi2 %.3f R2 %.4f'%(r['dip'],r['quad'],r['cold'],r['rms'],r['left'],r['leftu'],r['chi2'],r['R2']))
pickle.dump({k:{a:b for a,b in v.items() if a not in('x','af')} for k,v in R.items()},open('summary.pkl','wb'))
i=0; best=min([k for k in R if k[0]==i],key=lambda k:R[k]['rms']); print('best',best)
plt.rcParams.update({'font.size':13}); xx=np.arange(361)
fig,ax=plt.subplots(figsize=(11,6))
for k,c,ls,nm in [((i,0.0),'C0','-','LTB (no shell dipole)'),(best,'C3','--',f'Szekeres, shell dipole at {abs(best[1])*100:.0f}% of the no-crossing limit')]:
    ax.plot(xx,R[k]['x'],color=c,ls=ls,lw=2.2,label=f"{nm}: dipole {R[k]['dip']:.4f}, RMS {R[k]['rms']:.4f}, leftover {R[k]['left']:.4f}")
ax.plot(xx,full,'k',lw=1,label=f'observed (Planck SMICA + dipole), RMS {rms(full):.2e}')
ax.axhline(0,color='k',lw=.4); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_ylabel('ΔT/T'); ax.set_xlabel('direction on the ring (deg): 0 = away from CoE, 180 = toward CoE')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../szekeres_cmb_ltb_vs_szekeres.png',dpi=150)
fig,ax=plt.subplots(figsize=(11,6)); ax.errorbar(TH1[IMIG],Y,yerr=[MLO,MHI],fmt='o',ms=4,color='k',lw=.8,label='Migkas et al. 2020, A/A_all (71 bins)')
for k,c,ls,nm in [((i,0.0),'C0','-','LTB (no shell dipole)'),(best,'C3','--','Szekeres (same as CMB plot)')]:
    ax.plot(TH1,R[k]['af'],color=c,ls=ls,lw=2.2,label=f"{nm}: χ² {R[k]['chi2']:.2f}, R² {R[k]['R2']:.3f}")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel('Galactic longitude l (deg)'); ax.set_ylabel('A / A_all')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../szekeres_migkas_ltb_vs_szekeres.png',dpi=150)
