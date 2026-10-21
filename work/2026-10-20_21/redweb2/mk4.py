import json,numpy as np;exec(open('plot4.py').read())
T=json.load(open('thin4_results.json'));C=json.load(open('densweb/cmb_results.json'));sm=('smooth red (no web)',T['smooth']['Dl'],T['smooth']['chi2'])
cur=[sm];print('THIN4 f best mean r2best')
for f,L in T['thin'].items():
  c=[d['chi2'] for d in L];b=L[int(np.argmin(c))];cur.append((f'kept {f} (best of 5)',b['Dl'],b['chi2']));print(f,'%.4g %.4g %.2f'%(min(c),np.mean(c),r2(b['Dl'])))
for s in ('0.01','0.0175','0.02'):d=T['scaled_'+s];cur.append((f'strength {s} (all density-kept bodies)',d['Dl'],d['chi2']))
fig(cur,'thin4_log.png','Density-shaped thinning, web read at 4 Gpc')
F={};
for nm in ('orig','densweb_s1','densweb_s2','densweb_s3'):
  F[nm]={k:dict(chi2=v['chi2'],r2=r2(v['Dl']),Dl=v['Dl']) for k,v in C[nm+'_R4.0'].items()}
F['smooth']=dict(chi2=T['smooth']['chi2'],r2=r2(T['smooth']['Dl']))
json.dump(F,open('force4_results.json','w'))
for nm,t in [('orig','seed-2000 web, all bodies, strength sweep, 4 Gpc'),('densweb_s1','density-following web seed 1, strength sweep, 4 Gpc')]:
  fig([sm]+[(f'strength {k}',v['Dl'],v['chi2']) for k,v in F[nm].items()],f'force4_log_{nm}.png',t)
print('FORCE4 s orig dw1 dw2 dw3 (chi2/r2)')
for k in F['orig']:print(k,' '.join('%.4g/%.2f'%(F[n][k]['chi2'],F[n][k]['r2']) for n in ('orig','densweb_s1','densweb_s2','densweb_s3')))
D={};SS=['1.0','0.1','0.02','0.0175','0.01']
for R in ('4.0','0.9'):
  cur=[sm]
  for s in (1,2,3):
    for k in SS:v=C[f'densweb_s{s}_R{R}'][k];D[f's{s}_R{R}_{k}']=dict(chi2=v['chi2'],r2=r2(v['Dl']),Dl=v['Dl']);cur.append((f'seed{s} strength {k}',v['Dl'],v['chi2']))
  fig(cur,'densweb_log.png' if R=='4.0' else 'densweb_log_R0.9.png',f'Density-following web, read at {R} Gpc')
  print('DENSWEB R',R);[print(k,'%.4g/%.2f'%(D[f's{s}_R{R}_{k}']['chi2'],D[f's{s}_R{R}_{k}']['r2'])) for s in (1,2,3) for k in SS]
json.dump(D,open('densweb_results.json','w'))
