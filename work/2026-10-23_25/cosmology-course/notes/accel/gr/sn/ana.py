import pickle,numpy as np,json
from scipy.optimize import least_squares
from scipy.integrate import quad
from astropy.coordinates import SkyCoord; import astropy.units as un
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
def dl_lcdm(z,Om): return (1+z)*quad(lambda x:1/np.sqrt(Om*(1+x)**3+1-Om),0,z)[0]
def dl_eds(z): return 2*(1+z)*(1-1/np.sqrt(1+z))
def fitq(z,d,zmax=0.1):
    m=(z>=0.01)&(z<=zmax); z,d=z[m],d[m]
    f=lambda p: np.log(z/p[0]*(1+(1-p[1])*z/2-(1-p[1]-3*p[1]**2+p[2])*z*z/6))-np.log(d)
    r=least_squares(f,[0.3,0.5,1.0]); return r.x
# Pantheon+
import pandas as pd
P=pd.read_csv('pantheon.dat',sep=r'\s+'); C=np.loadtxt('pantheon_cov.cov',skiprows=1); N=int(np.loadtxt('pantheon_cov.cov',max_rows=1)); C=C.reshape(N,N)
sel=(P.IS_CALIBRATOR==0)&(P.zHD>0.01)
gal=SkyCoord(ra=P.RA.values*un.deg,dec=P.DEC.values*un.deg).galactic
def unit(l,b): l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
A=unit(15.166,30.0); nsn=unit(gal.l.deg,gal.b.deg); ang=np.degrees(np.arccos(np.clip(A@nsn,-1,1)))  # 0 = away from CoE
out={}
zg=np.linspace(0.005,2.0,400)
for i,nm in ((0,'flattest valley point (D 46.7 Gpc, core 0.093)'),(2,'best-fit valley point (D 55.6 Gpc, core 0.094)')):
    d=pickle.load(open(f'dl_{i}.pkl','rb')); th=[];DL=[];Q=[]
    for b,a in sorted(d['rays'].items()):
        look=180.0-b; z,dA=a[:,0],a[:,1]; u,k=np.unique(z,return_index=True); dl=(1+u)**2*dA[k]
        H0,q0,j0=fitq(u,dl); H0b,q0b,_=fitq(u,dl,0.05)
        th.append(look); Q.append((look,H0,q0,j0,q0b,u.max())); DL.append(np.interp(zg,u,dl,right=np.nan))
    o=np.argsort(th); th=np.array(th)[o]; DL=np.array(DL)[o]; Q=[Q[k] for k in o]
    avg=np.nanmean(DL,0); H0a,q0a,j0a=fitq(zg,avg)
    # per-SN model mu
    s=sel.values.copy(); z=P.zHD.values; mm=np.full(len(z),np.nan)
    for k in np.where(s)[0]:
        col=np.array([np.interp(z[k],zg,DL[j]) for j in range(len(th))]); 
        if np.isnan(col).any(): continue
        mm[k]=np.interp(ang[k],th,col)
    out[i]=dict(name=nm,th=th,DL=DL,Q=Q,avg=avg,qavg=(H0a,q0a,j0a),mm=mm,minRp=d['minRp'])
use=sel.values&np.isfinite(out[0]['mm'])&np.isfinite(out[2]['mm'])
idx=np.where(use)[0]; Ci=np.linalg.inv(C[np.ix_(idx,idx)]); mb=P.m_b_corr.values[idx]; zz=P.zHD.values[idx]; one=np.ones(len(idx))
def chi2(dl):
    r=mb-5*np.log10(dl); a=one@Ci@r/(one@Ci@one); rr=r-a; return float(rr@Ci@rr)
res={'N_SN':int(len(idx)),'N_dropped':int(sel.sum()-len(idx)),'zmax_used':float(zz.max())}
res['LCDM_Om0.3']=chi2(np.array([dl_lcdm(z,0.3) for z in zz])); res['EdS']=chi2(dl_eds(zz))
for i in (0,2):
    res[f'model_{i}_per_direction']=chi2(out[i]['mm'][idx]); res[f'model_{i}_ring_average']=chi2(np.interp(zz,zg,out[i]['avg']))
    res[f'q_{i}']=[dict(look_angle=q[0],H0_over_global=q[1]/0.3,q0=q[2],j0=q[3],q0_fit_to_0p05=q[4],zmax_reached=q[5]) for q in out[i]['Q']]
    res[f'q_{i}_ring_avg']=dict(H0_over_global=out[i]['qavg'][0]/0.3,q0=out[i]['qavg'][1],j0=out[i]['qavg'][2])
json.dump(res,open('sn_results.json','w'),indent=1); pickle.dump(out,open('sn_out.pkl','wb'))
for k,v in res.items(): print(k,v if not isinstance(v,list) else '\n'+'\n'.join(str({a:round(b,4) for a,b in x.items()}) for x in v))
# plots
plt.rcParams.update({'font.size':13}); mu=lambda d:5*np.log10(d)
fig,ax=plt.subplots(figsize=(11,6)); E=dl_eds(zg); m0=np.isfinite(out[0]['avg'])
def off(dm): w=(zg>0.01)&(zg<0.1)&np.isfinite(dm); return dm-dm[w].mean()
ax.plot(zg,off(mu(out[0]['avg'])-mu(E)),'C3',lw=2.4,label='model, flattest valley point (ring average)')
ax.plot(zg,off(mu(out[2]['avg'])-mu(E)),'C1',lw=2.0,ls='--',label='model, best-fit valley point (ring average)')
L_=np.array([dl_lcdm(z,0.3) for z in zg]); ax.plot(zg,off(mu(L_)-mu(E)),'C0',lw=2.2,label='flat ΛCDM, Ωm = 0.3')
r=mb-mu(dl_eds(zz)); bins=np.geomspace(0.01,2.3,16); bc=[];bm=[];be=[]
w0=(zz>0.01)&(zz<0.1); r=r-r[w0].mean()
for a,b in zip(bins[:-1],bins[1:]):
    k=(zz>=a)&(zz<b)
    if k.sum()>3: bc.append(np.median(zz[k])); bm.append(r[k].mean()); be.append(r[k].std()/np.sqrt(k.sum()))
ax.errorbar(bc,bm,yerr=be,fmt='ko',ms=6,label='Pantheon+ (binned, offset free)')
ax.axhline(0,color='k',lw=.5); ax.set_xlim(0,2.05); ax.set_xlabel('redshift z'); ax.set_ylabel('distance modulus minus Einstein–de Sitter (mag)')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../sn_mu_residual_vs_eds.png',dpi=150)
fig,ax=plt.subplots(figsize=(11,6))
for i,c,ls in ((0,'C3','-o'),(2,'C1','--s')):
    Q=out[i]['Q']; ax.plot([q[0] for q in Q],[q[2] for q in Q],ls,color=c,lw=2,ms=7,label=out[i]['name'])
ax.axhline(0.5,color='k',ls=':',label='Einstein–de Sitter, q0 = 0.5'); ax.axhline(-0.55,color='C0',ls=':',label='flat ΛCDM Ωm 0.3, q0 = −0.55')
ax.set_xlim(0,180); ax.set_xticks(range(0,181,30)); ax.set_xlabel('viewing direction (deg): 0 = away from CoE, 180 = toward CoE'); ax.set_ylabel('inferred deceleration q0 (fit 0.01 < z < 0.1)')
ax.legend(loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,fontsize=11,frameon=False); fig.tight_layout(); fig.savefig('../sn_q0_vs_direction.png',dpi=150)
