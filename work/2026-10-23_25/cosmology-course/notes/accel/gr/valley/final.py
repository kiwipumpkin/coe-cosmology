import sys,json,pickle,numpy as np,warnings; warnings.filterwarnings('ignore')
sys.path.insert(0,'..'); sys.path.insert(0,'/workspace/cosmology-course/notes/paper_repro/grfit')
from grcore import stats,norm
from core import ltb as LTBfit, shifted, chi, r2, TH1, IMIG, Y, MLO, MHI
import healpy as hp, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, toy
P=json.load(open('valley_points.json')); T=pickle.load(open('toy_0_1_2_3_5_6_7.pkl','rb'))
T0=2.7255; u=lambda l,b: np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(15.166,30.0); dv=u(264.021,48.253); beta=3362.08e-6/T0; Pp=dv-(dv@A)*A; Pp/=np.linalg.norm(Pp)
th=np.radians(np.arange(361)); n=np.outer(np.cos(th),A)+np.outer(np.sin(th),Pp)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True); M=hp.read_map('/workspace/planck/mask.fits',field=0,nest=True); ns=hp.get_nside(I)
cmb=np.empty(361); mf=np.empty(361)
for i,v in enumerate(n): pix=hp.query_disc(ns,v,np.radians(1.0),nest=True); cmb[i]=I[pix].mean(); mf[i]=1-M[pix].mean()
bad=(mf>0.5)|(np.abs(np.degrees(np.arcsin(n[:,2])))<10); full=cmb/T0+beta*(n@dv); full-=full[:-1].mean()
rms=lambda x,m=None: np.sqrt(np.mean((x[:-1] if m is None else x[:-1][~m[:-1]])**2))
rows=[]
for i in sorted(T):
    p=P[i]; r=T[i]; x=norm(r['x']); s=stats(x)
    L=toy.make(p)[0].func(*toy.make(p)[0].args[:3],r['rmin']); R=L.Rs(L.rg,L.tg,dx=1); tt=(L.tg>=L.t1)&(L.tg<=L.t1+r['dt']+0.5)
    rows.append(dict(i=i,**{k:p[k] for k in ('D','mcc','kin','kout','l0','Rin','Rout','chi2','R2','Dtoy')},view=r['dt'],rmin=r['rmin'],minRp_used=float(R[:,tt].min()),dip=s[0],quad=s[1],cold=float(x.min()),rms=rms(x),left=rms(x-full),left_unm=rms(x-full,bad),x=x.tolist()))
    print({k:(round(v,4) if isinstance(v,float) else v) for k,v in rows[-1].items() if k!='x'},flush=True)
json.dump(rows,open('valley_cmb.json','w'),indent=1)
rows.sort(key=lambda r:r['D']); best=min(rows,key=lambda r:r['chi2']); flat=min([r for r in rows if r['chi2']<min(q['chi2'] for q in P)+1],key=lambda r:r['rms'])
print('BEST',best['i'],'FLATTEST(dchi2<1)',flat['i'])
plt.rcParams.update({'font.size':13})
fig,ax=plt.subplots(figsize=(11,6)); g0=min(q['chi2'] for q in P)
for lab,sel,mk in [('Migkas χ² within 1 of best',lambda r:r['chi2']<g0+1,'o'),('Migkas χ² within 4 of best',lambda r:g0+1<=r['chi2']<g0+4,'s')]:
    q=[r for r in rows if sel(r)]; ax.plot([r['D'] for r in q],[r['rms'] for r in q],mk,ms=10,label=lab)
for r in rows: ax.annotate(f"core {r['mcc']:.3f}, χ² {r['chi2']:.1f}",(r['D'],r['rms']),textcoords='offset points',xytext=(6,6),fontsize=10)
ax.axhline(rms(full),ls='--',color='k',label=f'observed sky on this ring, RMS {rms(full):.2e}')
ax.set_xlabel('distance to the CoE, D (Gpc)'); ax.set_ylabel('calculated CMB ΔT/T RMS'); ax.set_ylim(0,None)
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../valley_cmb_rms.png',dpi=150)
fig,ax=plt.subplots(figsize=(11,6)); ax.errorbar(TH1[IMIG],Y,yerr=[MLO,MHI],fmt='o',ms=4,color='k',lw=.8,label='Migkas et al. 2020, A/A_all (71 bins)')
for r,c,nm in [(best,'C0','best Migkas fit'),(flat,'C3','flattest CMB in the χ²<best+1 valley')]:
    m=shifted(LTBfit(r['D'],r['mcc'],r['kin'],r['kout']),r['l0']); ax.plot(TH1,m,color=c,lw=2.2,ls='-' if c=='C0' else '--',label=f"{nm}: D {r['D']:.1f} Gpc, core {r['mcc']:.3f}, χ² {chi(m):.2f}, R² {r2(m):.3f}")
ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_xlabel('Galactic longitude l (deg)'); ax.set_ylabel('A / A_all')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../valley_migkas_best_vs_flat.png',dpi=150)
fig,ax=plt.subplots(figsize=(11,6)); xx=np.arange(361)
for r,c,nm in [(best,'C0','best Migkas fit'),(flat,'C3','flattest CMB in valley')]:
    ax.plot(xx,r['x'],color=c,lw=2.2,ls='-' if c=='C0' else '--',label=f"{nm}: RMS {r['rms']:.4f}, leftover vs Planck RMS {r['left']:.4f}")
ax.plot(xx,full,color='k',lw=1,label=f'observed (Planck SMICA + dipole), RMS {rms(full):.2e}')
ax.axhline(0,color='k',lw=.4); ax.set_xlim(0,360); ax.set_xticks(range(0,361,45)); ax.set_ylabel('ΔT/T')
ax.set_xlabel('direction on the ring (deg): 0 = away from CoE, 180 = toward CoE')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../valley_cmb_best_vs_flat.png',dpi=150)
