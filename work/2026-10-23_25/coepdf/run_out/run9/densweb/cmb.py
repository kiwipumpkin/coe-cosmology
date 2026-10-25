import sys,json,numpy as np,os
src=open('/workspace/redweb2/thin4_run.py').read()
head=src[:src.index('web=G.generate')];body=src[src.index('# full-sky generalisation'):src.index('res=dict(')]
g={};exec(head,g);exec(body,g)
D='/workspace/redweb2/densweb/';out=D+'cmb_results.json'
res=json.load(open(out)) if os.path.exists(out) else {}
res['smooth']=dict(Dl=g['Dl'](g['sm']),chi2=g['chi'](g['Dl'](g['sm'])))
S=[.01,.0175,.02,.05,.1,.2,.3,.4,.5,.6,.7,.8,.9,1.0]
webs=[('orig',D+'orig.npz')]+[(f'densweb_s{s}',D+f'densweb_s{s}.npz') for s in (1,2,3)]
for nm,f in webs:
  z=np.load(f);g['X'],g['M']=(z['arr_0'],z['arr_1']) if nm=='orig' else (z['X'],z['M'])
  for R in (4.0,0.9):
    key=f'{nm}_R{R}'
    if key in res:continue
    w=g['sky'](np.full(len(g['V']),R/g['RV']));np.save(D+key+'_w.npy',w)
    res[key]={str(s):dict(Dl=(d:=g['Dl'](g['sm']+s*w)),chi2=g['chi'](d)) for s in S}
    json.dump(res,open(out,'w'));print(key,{k:round(v['chi2'],1) for k,v in res[key].items()},flush=True)
