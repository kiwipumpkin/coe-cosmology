"""Shadow width: rays (model angle theta from anti-centre; CoE at 180) that end without reaching 3000 K with fixed K.
usage: shadow.py law mode vm ; env XOVR, KFIX, OUTTAG, LMSET (override lM)"""
import os,sys,json,numpy as np
src=open('cmb4.py').read()
law,mode,vm=sys.argv[1:4]; sys.argv=sys.argv[:4]+['30']
x=json.loads(os.environ['XOVR'])
if 'LMSET' in os.environ:
    import re
    names=['lM','u','p','lf','lL','s'] if mode=='fwd' else None
    i=(names.index('lM') if mode=='fwd' else (['lv','p','lM','lf','lL','s'] if vm=='gtc' else ['vo','p','lM','lf','lL','s']).index('lM')); x[i]=float(os.environ['LMSET']); os.environ['XOVR']=json.dumps(x)
os.environ['OUTD']='gr2/v4/shadow'
pre=src.split("t0=time.time(); ells=")[0]
exec(pre)
K=float(os.environ['KFIX']); lnth=3*np.log(3000/K)
def probe(e):
    rd=trace(np.radians(e)); q=recpoint(rd,lnth)
    return dict(ell=float(e),rec=q is not None,stop=rd['stop'],phi_e_deg=float(np.degrees(q['phi'])) if q else None,ncross=q['nc'] if q else None,
                T_K=float(3000*np.exp(-q['lnz1'])) if q else None,zmax=float(np.exp(np.nanmax(rd['lnz1']))-1))
grid=[probe(e) for e in np.arange(120,240.01,2.)]
def edge(a,b,ra):   # bisect between rec state ra at a and opposite at b
    for _ in range(5):
        m=0.5*(a+b); p=probe(m)
        if p['rec']==ra: a=m
        else: b=m
    return 0.5*(a+b)
edges=[]
for g1,g2 in zip(grid[:-1],grid[1:]):
    if g1['rec']!=g2['rec']: edges.append(edge(g1['ell'],g2['ell'],g1['rec']))
nshadow=[g['ell'] for g in grid if not g['rec']]
width=(edges[-1]-edges[0]) if len(edges)>=2 else (0.0 if not nshadow else (2.0 if len(edges)<2 else None))
# near-edge rays: deflection of emitter azimuth and T
near=[g for g in grid if g['rec'] and nshadow and min(abs(g['ell']-s) for s in nshadow)<=6]
out=dict(tag=OT,law=law,mode=mode,vm=vm,x=x,names=nm,M=float(P.get('M',np.nan)),Uo_c=float(Uo),D_Gpc=float(Ro*lam),chi2=chi,l0=l0,K=K,edges=edges,shadow_width_deg=width,
  shadow_dirs_2deg=nshadow,n_shadow_2deg=len(nshadow),grid=grid,near_edge=near)
json.dump(out,open(f'gr2/v4/shadow/{OT}.json','w'),indent=1,default=lambda o: np.asarray(o).tolist()); print('DONE',OT,'M %.3g'%out['M'],'width',width,'edges',edges,'chi2 %.2f'%chi,flush=True)
