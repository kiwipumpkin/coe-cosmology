# Inverse fit with rings: delta(l,b)=sum_i A_i exp(-0.5((theta_n(l,b)-theta_i)/w_i)^2), theta_n = angle from axis n.
# Target as A_required_web_solve.py: dipole-removed share of connect-the-dots y+strip_avg(delta) (71 Migkas bins,
# |dl|<45, all b, cos b) = Planck TT (2l+1)C_l share for m=l=2..35, total m=2..35 power of y fixed. Small norm penalty lam on delta.
# Grid 3 deg, FFT 4096 points (m<=35 targeted). Cases: shared axis fixed at CoE (180,-15), fixed at AoE (260,60), shared fitted axis; free axis per ring. N=2,4,6,8. Multi-start.
import numpy as np,json
from scipy.optimize import least_squares
rng=np.random.default_rng(0)
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
o=np.argsort(d['l_deg']);l=d['l_deg'][o];y=d['A_over_Aall'][o]
Ng=2**12;g=np.arange(Ng)*360/Ng
def spec(v):return 4*np.abs(np.fft.rfft(np.interp(g,l,v,period=360))/Ng)**2
p=np.loadtxt('COM_PowerSpect_CMB-TT-full_R3.01.txt');L=p[:,0].astype(int);Pl=(2*L+1)*2*np.pi*p[:,1]/(L*(L+1))
M=35;pl=Pl[(L>=2)&(L<=M)];pl/=pl.sum();T=spec(y)[2:M+1].sum()
lg=np.radians(np.arange(1.5,360,3));bg=np.radians(np.arange(-88.5,90,3));LL,BB=np.meshgrid(lg,bg);W=np.cos(BB)
xyz=np.stack([np.cos(BB)*np.cos(LL),np.cos(BB)*np.sin(LL),np.sin(BB)]).reshape(3,-1)
St=np.array([((np.abs(np.angle(np.exp(1j*(LL-np.radians(c)))))<np.radians(45))*W).ravel() for c in l]);St/=St.sum(1,keepdims=True)
def vec(lo,la):lo,la=np.radians(lo),np.radians(la);return np.array([np.cos(la)*np.cos(lo),np.cos(la)*np.sin(lo),np.sin(la)])
def sky(rings):
  s=np.zeros(xyz.shape[1])
  for (al,ab,th,w,A) in rings:
    t=np.degrees(np.arccos(np.clip(vec(al,ab)@xyz,-1,1)));s+=A*np.exp(-0.5*((t-th)/w)**2)
  return s
lam=1e-3
def unpack(x,N,mode,ax):
  if mode=='free':r=x.reshape(N,5);return [tuple(q) for q in r]
  a=ax if ax is not None else x[:2];q=(x if ax is not None else x[2:]).reshape(N,3);return [(a[0],a[1],*qq) for qq in q]
def res(x,N,mode,ax):
  s=sky(unpack(x,N,mode,ax));D=St@s;P=spec(y+D)[2:M+1]
  return np.concatenate([(P-T*pl)/T,[lam*s.std()]])
def fit(N,mode,ax,starts=4):
  best=None
  for k in range(starts):
    rp=[[rng.uniform(5,175),rng.uniform(2,30),rng.normal(0,.5)] for _ in range(N)]
    if mode=='free':x0=np.array([[rng.uniform(0,360),np.degrees(np.arcsin(rng.uniform(-1,1)))]+r for r in rp]).ravel();lb=np.tile([-np.inf,-90,0,.5,-np.inf],N);ub=np.tile([np.inf,90,180,90,np.inf],N)
    else:
      x0=np.array(sum(rp,[]));lb=np.tile([0,.5,-np.inf],N);ub=np.tile([180,90,np.inf],N)
      if ax is None:x0=np.r_[rng.uniform(0,360),np.degrees(np.arcsin(rng.uniform(-1,1))),x0];lb=np.r_[-np.inf,-90,lb];ub=np.r_[np.inf,90,ub]
    r=least_squares(res,x0,bounds=(lb,ub),args=(N,mode,ax),max_nfev=60*len(x0),x_scale='jac')
    if best is None or r.cost<best.cost:best=r
  rings=unpack(best.x,N,mode,ax);s=sky(rings);D=St@s;P=spec(y+D);a=P[2:M+1]/P[2:M+1].sum()
  return dict(rings=[[float(v) for v in q] for q in rings],resid=float(np.sqrt(((a-pl)**2).mean())),drms=float(s.std()),dpeak=float(np.abs(s).max()),Drms=float(D.std()))
out={}
for name,mode,ax in [('CoE','shared',(180,-15)),('AoE','shared',(260,60)),('shared_fit','shared',None),('free','free',None)]:
  for N in [2,4,6,8]:
    r=fit(N,mode,ax);out[f'{name}_{N}']=r;print(name,N,'resid %.4f drms %.3f peak %.3f Drms %.3f'%(r['resid'],r['drms'],r['dpeak'],r['Drms']),flush=True)
json.dump(out,open('A_required_rings_results.json','w'),indent=1)
