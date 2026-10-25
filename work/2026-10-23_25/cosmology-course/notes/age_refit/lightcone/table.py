import json,glob,numpy as np
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
rows=[]
for f in sorted(glob.glob('res_*.json')):
  r=json.load(open(f)); A=r['cmb_A']; ok=[q for q in A.get('rows',[]) if 't_e' in q]
  st=A.get('status_counts',{}); p=', '.join(f'{n}={v:.3g}' for n,v in zip(r['names'],r['x']))
  rows.append([r['law'],p,f"{r['M_kg']:.2g}",f"{r['D_Gpc']:.2f}",f"{r['U_c']:.2f}",f"{r['Hpar']:.0f}/{r['Hper']:.0f}",f"{r['age_Gyr']:.1f}",f"{r['chi2']:.1f}",f"{r['sigma']:.1f}",f"{r['l0']:.0f}",
   f"{A['K']:.2f}" if A.get('K') else 'none',str([round(v,2) for v in A['T_range']]) if 'T_range' in A else '-',f"{A['dTT_rms_m1rem']:.3f}" if 'dTT_rms_m1rem' in A else '-',
   f"{A['resid_rms']:.3f}" if 'resid_rms' in A else '-',f"{A['cancel_pct']:.0f}%" if 'cancel_pct' in A else '-',f"{A['corr_with_Gm1']:+.2f}" if 'corr_with_Gm1' in A else '-',
   f"{A.get('n_emit',0)}/60; {st}; beyondCoE {A.get('beyond_CoE_emitters','-')}",str(A.get('order',{}).get('zero_first','-'))])
hdr=['law','params','M kg','D Gpc','U/c','H par/perp','age','chi2','sigma','l0','K needed','T range K','dT/T rms','resid','cancel','corr','directions','0 first']
fig,ax=plt.subplots(figsize=(30,3+1.2*len(rows)));ax.axis('off')
t=ax.table(cellText=rows,colLabels=hdr,loc='center');t.auto_set_font_size(False);t.set_fontsize(11);t.scale(1,3);t.auto_set_column_width(range(len(hdr)))
fig.savefig('summary_table.png',dpi=80,bbox_inches='tight'); json.dump(dict(hdr=hdr,rows=rows),open('summary_table.json','w'),indent=1)
for r in rows: print(r)
