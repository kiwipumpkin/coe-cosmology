import json,glob
for f in sorted(glob.glob('res_*.json')):
  r=json.load(open(f)); print(f,{k:(round(v,4) if isinstance(v,float) else v) for k,v in r.items() if k in ['names','x','chi2','sigma','l0','U_c','Hpar','Hper','D_Gpc','age_Gyr','M_kg','rs_over_D']})
  from collections import Counter; print(' rays',Counter(v['stop'] for v in r['rays'].values()))
  for tm in 'AB':
    d=r['cmb_'+tm]; print(' ',tm,{k:d.get(k) for k in ['K','tstar_Gyr','n_emit','status_counts','beyond_CoE_emitters','T_range','dTT_rms_m1rem','target_rms_m1rem','resid_rms','cancel_pct','corr_with_Gm1','order']})
    ok=[q for q in d.get('rows',[]) if 't_e' in q]
    if ok: print('   z range',min(2.718**q['lnz1']-1 for q in ok),max(2.718**q['lnz1']-1 for q in ok),'Temit',min(q['T_emit'] for q in ok),max(q['T_emit'] for q in ok))
    pv=d.get('peculiar_velocity_variants')
    if pv: print('   pv',{b:(round(v['resid_rms'],3),round(v['cancel_pct']),round(v['corr_with_Gm1'],2)) for b,v in pv.items()})
