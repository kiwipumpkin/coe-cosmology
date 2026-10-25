import json,glob,os,numpy as np,importlib.util as u
s=u.spec_from_file_location('cc','/workspace/clump_scan/cmb-cancel.py');m=u.module_from_spec(s);s.loader.exec_module(m)
rows=[]
for f in sorted(glob.glob('gr2/v4/mscan/lM*.json')):
    d=json.load(open(f)); t=d['outtag']; sp=d.get('score_partial')
    if not sp: rows.append((t,'no score',d.get('K_note'))); continue
    gl=np.array(sp['gal_l']); g=np.interp(gl,np.arange(361),np.append(m.G_FULL,m.G_FULL[0]))-1
    sh=f'gr2/v4/mscan/shs_{t}.json'; W=pk=None
    if os.path.exists(sh):
        S=json.load(open(sh)); W=S['shadow_width_deg']; pk=max([x['T_K'] for x in S['grid'] if x['rec']])
    nr=[x['ell'] for x in d['rows'] if x['status']!='rec']
    rows.append(dict(tag=t,M=d['M'],chi2=d['chi2'],sigma=d['sigma'],Uo=d['Uo_c'],D=d['D_Gpc'],K=d['K_local_today_K'],cov=d['coverage'],shadow=W,shadow6=len(nr)*6,ringpk=pk,Tmax6=d['T_obs_range_K'][1],dTT=sp['model_dTT_rms_m1rem'],corr=float(np.corrcoef(sp['dTT'],g)[0,1]),res=sp['res_1e5']*1e-5,canc=sp['cancel'],ord=[x.get('t_e_Gyr') for x in d['rows'] if x['ell'] in (0.0,90.0,150.0)]))
json.dump(rows,open('gr2/v4/mscan/table.json','w'),indent=1,default=float)
for r in rows: print(r if isinstance(r,tuple) else ' '.join(f'{k}={v:.4g}' if isinstance(v,float) else f'{k}={v}' for k,v in r.items()))
