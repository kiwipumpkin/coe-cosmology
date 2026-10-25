import json,numpy as np,os,sys
os.environ['VMAX']='1'; sys.argv=['x','eta']
src=open('fit_forward.py').read().split('x=fit()')[0]; exec(src)
z=json.load(open('fwd/eta_vc_ckpt_free.json'))['x']; e=ev(np.array(z),True); o=e['o']; L=e['L']; lam=L['Hbar']/H70
D=L['D']; out=dict(lam=lam,D=D*lam,vo_c=L['Uo']/F.C,rows=[])
for f in (0.25,0.5,0.75,0.9):
    R=f*D; c=F.cmb_rrec(o,R); g=[F.gap_dr(o,R,d*D) for d in (0.01,0.05,0.1)]
    dtc,dto=F.dt_ls_formula(D*lam,L['Uo'],R*lam)
    th=np.arange(360); h=c['hit']
    out['rows'].append(dict(Rrec_over_D=f,Rrec_Gpc=R*lam,psi_crit_deg=c['psi_crit_deg'],n_dirs_blocked=int(h.sum()),
      R_critical_Gpc=c['R_critical_dist_from_us']*lam,t_obs_cross_Gyr=c['t_obs_cross']*lam,
      tLS_180=c['tLS'][180]*lam,tLS_170=c['tLS'][170]*lam,tLS_edge=float(np.nanmin(c['tLS']))*lam,
      gap180_int=c['gap_180_integrated']*lam,gap180_constU=c['gap_formula_constU']*lam,gap180_intformula=c['gap_formula_int']*lam,
      gap_dr=[(d,a*lam,b*lam) for d,(a,b) in zip((0.01,0.05,0.1),g)],dt_formula_CoE=dtc,dt_formula_obs=dto,
      Dm_180=float(c['Dm'][180]),tLS_theta=(c['tLS']*lam).tolist()))
json.dump(out,open('fwd/eta_vc_rrec.json','w'),indent=1,default=float)
for r in out['rows']: print({k:v for k,v in r.items() if k!='tLS_theta'})
