import pickle,json,numpy as np,healpy as hp
d=pickle.load(open('bf.pkl','rb'));R,p=d['R'],d['params'];AR=0.0045;TH=np.arange(0,360,2.0);T0=2.7255
u=lambda l,b:np.array([np.cos(np.radians(b))*np.cos(np.radians(l)),np.cos(np.radians(b))*np.sin(np.radians(l)),np.sin(np.radians(b))])
A=u(p['l0'],30.);dv=u(264.021,48.253);Pp=dv-(dv@A)*A;Pp/=np.linalg.norm(Pp);th=np.radians(TH);n=np.outer(np.cos(th),A)+np.outer(np.sin(th),Pp)
I=hp.read_map('/workspace/planck/smica.fits',field=0,nest=True);ns=hp.get_nside(I)
ob=np.array([I[hp.query_disc(ns,v,np.radians(1.),nest=True)].mean() for v in n])/T0+3362.08e-6/T0*(n@dv);ob-=ob.mean()
w=R[p['f']]['res'][AR];o=R[0.0]['res'][AR];ok=np.isfinite(w[:,0])&np.isfinite(o[:,0])
S=lambda x:dict(mean=float(np.mean(x)),mean_abs=float(np.mean(abs(x))),median_abs=float(np.median(abs(x))),rms=float(np.sqrt(np.mean(x**2))),range=float(x.max()-x.min()))
dw=w[ok,0]/w[ok,0].mean()-1;do=o[ok,0]/o[ok,0].mean()-1
out=dict(a_rec=AR,one_plus_z_approx=1/AR,missing=int((~ok).sum()),missing_angles=TH[~ok].tolist(),with_mass=S(dw),no_mass=S(do),observed=S(ob),
 Tobs_mean_K_with=float(3000*w[ok,0].mean()),Tobs_mean_K_no=float(3000*o[ok,0].mean()),R_rec_Gpc_with=[float(np.nanmin(w[:,2])*float(np.exp(json.load(open('fit_soft.json'))['best']['ln_D_Gpc']))/p['uD']),float(np.nanmax(w[:,2])*float(np.exp(json.load(open('fit_soft.json'))['best']['ln_D_Gpc']))/p['uD'])],params=p)
for k in('rms','range','mean_abs','median_abs'): out['pct_change_adding_mass_'+k]=100*(out['with_mass'][k]/out['no_mass'][k]-1);out['pct_change_removing_mass_'+k]=100*(out['no_mass'][k]/out['with_mass'][k]-1)
json.dump(out,open('density_cmb.json','w'),indent=1,default=float);print(json.dumps({k:v for k,v in out.items() if k!='params'},indent=1,default=float))
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt;plt.rcParams.update({'font.size':15})
fig,ax=plt.subplots(figsize=(13,8));lab=lambda s:f"RMS={s['rms']:.2e}, range={s['range']:.2e}, mean|·|={s['mean_abs']:.2e}, median|·|={s['median_abs']:.2e}"
for y,nm,s,c in[(dw,'with mass (soft core+tail, 3% lock)',out['with_mass'],'C0'),(do,'mass removed',out['no_mass'],'C3')]:
    Y=np.full(len(TH),np.nan);Y[ok]=y;ax.plot(TH,Y,c,lw=2.5,label=f'{nm}: '+lab(s))
ax.plot(TH,ob,'k',lw=2,label='Planck SMICA+dipole ring: '+lab(out['observed']))
ax.set_xlabel('ring angle (deg; CoE at 180°)');ax.set_ylabel('ΔT/T');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.0),frameon=False,fontsize=12)
fig.tight_layout();fig.savefig('density_cmb.png',dpi=110)
