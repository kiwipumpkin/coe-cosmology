import pickle,json,numpy as np
D=pickle.load(open('temp.pkl','rb'));obs=np.load('obs_ring.npy')[:-1][::2];TH=np.arange(0,360,2.0)
st=lambda d:dict(rms=float(np.sqrt(np.mean(d**2))),range=float(d.max()-d.min()),mean_abs=float(np.mean(abs(d))),med_abs=float(np.median(abs(d))))
ob=st(obs);rows=[];sel={}
for u in (1.05,3.0,10.5):
    m,n=D[(u,0.0178)],D[(u,0.0)];RV=46.7/u
    for ar in m['x']:
        opz=lambda r:1/(r['x'][ar]*ar)
        a,b=opz(m),opz(n);ok=np.isfinite(a)&np.isfinite(b)
        if ok.sum()<0.9*len(TH): continue
        da=a[ok]/a[ok].mean()-1;db=b[ok]/b[ok].mean()-1;sa,sb=st(da),st(db)
        res=float(np.sqrt(np.mean((da-obs[ok]+np.mean(obs[ok]))**2)))
        dt=(m['t1']-np.nanmean(m['te'][ar]))*RV*3.2616
        r=dict(speed=u,a_rec=ar,missing=int((~ok).sum()),with_mass=sa,no_mass=sb,
          drop_rms_pct=100*(1-sa['rms']/sb['rms']),drop_range_pct=100*(1-sa['range']/sb['range']),
          model_over_obs_range=sa['range']/ob['range'],resid_rms_vs_planck=res,
          time_since_recomb_Gyr=dt,R_rec_Gpc=[float(np.nanmin(m['Rem'][ar])*RV),float(np.nanmax(m['Rem'][ar])*RV)],
          T_today_local_K=3000*ar,mean_Tobs_K=float(3000*np.nanmean(a)))
        rows.append(r)
    best=min([r for r in rows if r['speed']==u],key=lambda r:r['resid_rms_vs_planck']);sel[u]=best
json.dump(dict(obs=ob,rows=rows,best=sel),open('cmb_results.json','w'),indent=1)
print('obs',ob)
for r in rows: print(r['speed'],r['a_rec'],r['missing'],'rmsM %.4g rmsN %.4g drop %.0f%% / %.0f%% resid %.4g dt %.3g Gyr R %s'%(r['with_mass']['rms'],r['no_mass']['rms'],r['drop_rms_pct'],r['drop_range_pct'],r['resid_rms_vs_planck'],r['time_since_recomb_Gyr'],np.round(r['R_rec_Gpc'],2)))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':15})
for u,r in sel.items():
    m,n=D[(u,0.0178)],D[(u,0.0)];ar=r['a_rec'];a,b=1/(m['x'][ar]*ar),1/(n['x'][ar]*ar);ok=np.isfinite(a)&np.isfinite(b)
    fig,ax=plt.subplots(figsize=(13,8))
    for y,lab,s,c in[(a/np.nanmean(a[ok])-1,'with fitted core mass',r['with_mass'],'C0'),(b/np.nanmean(b[ok])-1,'no mass',r['no_mass'],'C3'),(obs,'Planck SMICA+dipole ring',ob,'k')]:
        ax.plot(TH,np.where(ok|(c=='k'),y,np.nan),c,lw=2.5,label=f"{lab}: RMS={s['rms']:.2e}, range={s['range']:.2e}")
    ax.set_xlabel('ring angle (deg; CoE at 180°)');ax.set_ylabel('ΔT/T');ax.set_title(f'D=46.7 Gpc, space speed {u}c, recombination scale {ar} (best vs Planck); RMS drop with mass {r["drop_rms_pct"]:.0f}%',fontsize=12)
    ax.legend(loc='lower center',bbox_to_anchor=(.5,1.05),frameon=False);fig.tight_layout();fig.savefig(f'cmb_point_{u}c.png',dpi=110)
