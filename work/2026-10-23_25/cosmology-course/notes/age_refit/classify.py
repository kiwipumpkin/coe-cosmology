import json,glob,sys,numpy as np
files=sys.argv[1:] or sorted(glob.glob('gr2/rec3000_*.json'))
allc={}
for f in files:
    d=json.load(open(f)); rows=d['rows']; cls=[]
    for r in rows:
        if r['status']=='rec':
            ph=abs(r['phi_e_deg']); order=1+int(ph//180)+r['ncross']
            if r['ncross']>0 or ph>90: c='b_lensed_behind_core'
            elif r['R_e_over_D']<1 and 90<r['ell']<270: c='a_core_surface'
            else: c='outer_surface'
        else: c='c_no_solution:'+('trapped_W2<0' if r.get('diag',{}).get('W2',1)<=0 else ('A<=0_shell_crossing' if r.get('diag',{}).get('A',1)<=0 else r['status']))
        r['class']=c; r['image_order']=order if r['status']=='rec' else None; cls.append(c)
    from collections import Counter; cnt=Counter(cls)
    ok=[r for r in rows if r['status']=='rec']; s=dict(tag=d['tag'],counts=dict(cnt),Uo_global_c=d['Uo_c'],U_rel_local_flow_c=0.0)
    if ok:
        z=np.array([r['z'] for r in ok]); s['T_obs_K']=(float((3000/(1+z)).min()),float((3000/(1+z)).max()))
        a=[r for r in ok if r['class']=='a_core_surface']; o=[r for r in ok if r['class']=='outer_surface']
        if a and o: s['core_vs_outer_T_ratio']=float(np.mean([3000/(1+r['z']) for r in a])/np.mean([3000/(1+r['z']) for r in o]))
        s['t_e_Gyr']=(min(r['t_e_Gyr'] for r in ok),max(r['t_e_Gyr'] for r in ok)); s['flight_Gyr']=(min(r['flight_Gyr'] for r in ok),max(r['flight_Gyr'] for r in ok))
        s['near180']=[(r['ell'],r['class'],round(3000/(1+r['z']),5) if r['status']=='rec' else None) for r in rows if 150<=r['ell']<=210]
    for k in ('raw','nodip'):
        if k in d: s[k]=dict(cancel=round(d[k]['cancel_l0'],1),res_1e5=round(d[k]['res_l0_1e5'],3),m1=d[k]['modes'].get('1'),m2=d[k]['modes'].get('2'))
    s['chi2']=d['chi2']; s['sigma']=d['sigma']; allc[d['tag']]=s; print(json.dumps(s,default=float))
    json.dump(d,open(f,'w'),indent=1,default=float)
json.dump(allc,open('gr2/classify_summary.json','w'),indent=1,default=float)
