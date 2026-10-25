import os,sys,json,numpy as np
src=open('cmb4.py').read(); sys.argv=sys.argv[:4]+['30']; os.environ['OUTD']='gr2/v4/m253'
exec(src.split("t0=time.time(); ells=")[0])
K=float(os.environ['KFIX']); lnth=3*np.log(3000/K); res=[]
for e in np.arange(float(os.environ.get('RLO',160)),float(os.environ.get('RHI',200))+0.01,0.5):
    rd=trace(np.radians(e)); q=recpoint(rd,lnth)
    res.append(dict(ell=float(e),rec=q is not None,T_K=float(3000*np.exp(-q['lnz1'])) if q else None,phi_e_deg=float(np.degrees(q['phi'])) if q else None,t_e_Gyr=float((q['t']-tbo)/F.C*lam) if q else None))
json.dump(dict(K=K,l0=l0,rows=res),open(f'gr2/v4/m253/{OT}.json','w'),indent=1); print('DONE')
