import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.signal import find_peaks,savgol_filter
import half_fit as H;import acoustic_fit as A;import acoustic_model as M
ll=np.arange(M.LMAX+1);Lt=A.Lt;Pf=A.Pf;Pb=A.Pb;J=H.J
def peaks(D,lo=100):
    s=savgol_filter(D[lo:],61,3);pk,_=find_peaks(s,prominence=20);pk=(pk+lo)[:5];return pk.tolist(),[round(float(s[p-lo]),0) for p in pk]
def real(Df,seed=42):return Df*np.random.default_rng(seed).chisquare(2*np.maximum(ll,1)+1)/(2*np.maximum(ll,1)+1)
out={};models={}
Df=M.Dfield(J['free_W']['best']['p'])*H.CORR0;models['old']=(Df+M.Dfix,Df)
modes=[m for m in ('half','r2free') if 'best' in json.load(open(f'half_ckpt_{m}.json'))]
for m in modes:
    C=json.load(open(f'half_ckpt_{m}.json'));H.CORR[:]=C['corr'];Df=M.Dfield(C['best']['p'])*H.CORR;models[m]=(Df+M.Dfix,Df);out[m]=C['best']
bm=min(modes,key=lambda m:out[m]['l2_2500']['chi2_planck'])
sPl=savgol_filter(Pf[:,1][98:],61,3);pk,_=find_peaks(sPl,prominence=20);pk=pk[:5];out['peaks_planck']=((pk+100).tolist(),[round(float(sPl[p]),0) for p in pk])
for k,(D,Df) in models.items():out.setdefault(k,{})['peaks']=peaks(D)
out['peaks_lcdm']=peaks(Lt);out['lcdm']=J['lcdm']
D,Df=models[bm];Dr=M.Dfix+real(Df);st={k:A.stats(Dr,Df,*r) for k,r in dict(l2_30=(2,30),l2_2500=(2,2500)).items()};out[bm]['realization_seed42']=st
fig,ax=plt.subplots(2,1,figsize=(12,10))
for a,(l0,l1),key in zip(ax,((2,30),(2,M.LMAX)),('l2_30','l2_2500')):
    m=(ll>=l0)&(ll<=l1);a.plot(ll[m],Lt[m],'k-',lw=1.8,label=f'ΛCDM  χ²(Planck err)={J["lcdm"][key]["chi2_planck"]:.1f}, Pearson r²={J["lcdm"][key]["r2"]:.3f}')
    cv=np.sqrt(2/(2*np.maximum(ll,1)+1))*Df;a.fill_between(ll[m],(D-cv)[m],(D+cv)[m],color='r',alpha=.2,label='±1σ cosmic variance')
    a.plot(ll[m],D[m],'r-',lw=.7,alpha=.35,label=f'CoE mean (faint)  χ²(Planck err)={out[bm][key]["chi2_planck"]:.1f}, Pearson r²={out[bm][key]["r2"]:.3f}')
    a.plot(ll[m],Dr[m],'r-',lw=1.0 if l1>30 else 1.6,marker='' if l1>30 else 'o',ms=3,label=f'CoE + damped half-scale term, one realization  χ²(Planck err)={st[key]["chi2_planck"]:.1f}, Pearson r²={st[key]["r2"]:.3f}')
    mm=(Pf[:,0]>=l0)&(Pf[:,0]<=l1)
    if l1==30:a.errorbar(Pf[mm,0],Pf[mm,1],yerr=[Pf[mm,2],Pf[mm,3]],fmt='o',color='grey',ms=4,label='Planck 2018');a.set_ylim(0,2500)
    else:a.plot(Pf[mm,0],Pf[mm,1],'.',color='grey',ms=1.5,alpha=.5,label='Planck 2018 (unbinned)');a.set_xscale('log');a.set_ylim(0,7000)
    a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.legend(loc='lower center',fontsize=7,framealpha=.85)
ax[0].set_title(f'CoE acoustic field + damped half-scale term ({bm}), R=4 Gpc');fig.tight_layout();fig.savefig('acoustic_half_cl.png',dpi=100);plt.close()
fig,a=plt.subplots(figsize=(12,4.8))
for k,c,lab in (('old','tab:orange','old CoE (free_W)'),(bm,'r','new CoE (+half-scale term)')):
    a.plot(Pb[:,0],(Pb[:,1]-A.binned(models[k][0]))/Pb[:,2],'o-',color=c,ms=3,label=f'Planck − {lab}')
a.plot(Pb[:,0],(Pb[:,1]-A.binned(Lt))/Pb[:,2],'ko-',ms=3,label='Planck − ΛCDM');a.axhline(0,color='grey');a.set_xlabel('ℓ (binned, linear)');a.set_ylabel('residual / σ_Planck');a.legend(loc='lower center',fontsize=8);fig.tight_layout();fig.savefig('acoustic_half_resid.png',dpi=100);plt.close()
# leftover period: Lomb-Scargle on new residual
from scipy.signal import lombscargle
r=(Pb[:,1]-A.binned(models[bm][0]))/Pb[:,2];per=np.linspace(100,1500,1400);pw=lombscargle(Pb[:,0],r-r.mean(),2*np.pi/per);out['leftover_period']=float(per[pw.argmax()]);out['resid_rms']={k:float(np.sqrt(np.mean(((Pb[:,1]-A.binned(v[0]))/Pb[:,2])**2))) for k,v in models.items()}
out['resid_rms']['lcdm']=float(np.sqrt(np.mean(((Pb[:,1]-A.binned(Lt))/Pb[:,2])**2)));out['best_mode']=bm
np.save('resid_half.npy',np.c_[Pb[:,0],r]);json.dump(out,open('acoustic_half.json','w'),indent=1)
print(json.dumps({k:(v if k not in modes else {kk:v[kk] for kk in ('p','l2_2500','l2_30','peaks')}) for k,v in out.items() if k!='lcdm'},indent=0,default=str))
