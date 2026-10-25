import json,glob,numpy as np,sys
rows=[]
for f in sorted(glob.glob(sys.argv[1] if len(sys.argv)>1 else 'gr2/v4/s2/*.json')):
    d=json.load(open(f))
    if 'failed' in d: rows.append((d.get('outtag',f),'FAIL',d['failed'])); continue
    sp=d.get('score_partial',{}); rows.append((d['outtag'],round(d['chi2'],2),round(d['sigma'],2),'K=%.3g'%d.get('K_local_today_K',np.nan),d['K_note'][:30],'cov %s/%s'%(d.get('coverage'),d['n_dir']),'Tmean %.3g'%d.get('T_obs_mean_K',np.nan),'T',np.round(d.get('T_obs_range_K',[np.nan,np.nan]),3).tolist(),'dTT_rms %.3g'%sp.get('model_dTT_rms_m1rem',np.nan),'res %.3g'%(sp.get('res_1e5',np.nan)*1e-5),'canc %.1f'%sp.get('cancel',np.nan),'ord0<180',d.get('order_0_before_180'),'Uo %.3g'%d['Uo_c'],'M %.3g'%d['M']))
for r in rows: print(*r)
