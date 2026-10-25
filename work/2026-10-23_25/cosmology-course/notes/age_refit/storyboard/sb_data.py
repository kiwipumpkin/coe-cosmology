import os,sys,json,numpy as np
os.chdir('/workspace/cosmology-course/notes/age_refit')
sys.path.insert(0,'/workspace/cosmology-course/notes/age_refit')
b=json.load(open('gr2/v4/best/best_eta_fwd_gtc_1.json'))
os.environ['XOVR']=json.dumps(b['x']); sys.argv=['cmb4.py','eta','fwd','gtc','30']
exec(open('cmb4.py').read().split("t0=time.time(); ells=")[0])
K=b['K_local_today_K']; lnth=3*np.log(3000/K)
xs=np.linspace(np.log(1e-3)+1,np.log(3e4)-1,1500)
tbmin=TB_MIN; ts=np.concatenate([tbmin+(tn-tbmin)*np.logspace(-6,0,40)])
Rg=np.array([L.R(t,xs) for t in ts]); rg=np.array([L.rho_rel(t,xs) for t in ts])
Hobs=np.array([L.Rdot(t,X)[0]/L.R(t,X)[0] for t in ts]); Robs=np.array([L.R(t,X)[0] for t in ts])
from multiprocessing import Pool
angs=[0,20,40,60,80,100,120,135,144,146,148,148.5,150,165,180]
def one(e):
    rd=trace(np.radians(e)); q=recpoint(rd,lnth); return e,rd,q
with Pool(15) as p: res=p.map(one,angs)
rays={}
for e,rd,q in res:
    rays[str(e)]=dict(t=rd['t'].tolist(),R=[float(L.R(tt,np.array([xx]))[0]) for tt,xx in zip(rd['t'],rd['x'])],phi=rd['phi'].tolist(),stop=rd['stop'],
      rec=None if q is None else dict(t=float(q['t']),R=float(L.R(q['t'],np.array([q['x']]))[0]),phi=float(q['phi']),z=float(np.exp(q['lnz1'])-1)))
np.savez('storyboard/sb_data.npz',ts=ts,xs=xs,R=Rg,rho=rg,H=Hobs,Robs=Robs,tn=tn,tbo=tbo,lam=lam,Ro=Ro,Uo=Uo,rho_o=rho_o,K=K,C=F.C,M=P['M'],
  mR=np.array([L.m(x_) for x_ in xs]))
json.dump(rays,open('storyboard/sb_rays.json','w'))
print('ok',lam,F.C,tn,tbo,Ro,Hobs[-1])
