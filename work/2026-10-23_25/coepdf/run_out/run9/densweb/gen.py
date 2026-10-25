# density-following random web: same bodies/mass spectrum as seed-2000 web (void -D/+W pairs move together, clusters alone),
# positions redrawn: uniform in observer ball 0.5<r<20 Gpc, accepted with prob shp(|P-PC|/42)/max  -> number (and mean mass) density ∝ rho(r_CoE)
import os,sys,numpy as np
A='/workspace/coe-cosmology-gh/draft/current/analyses';sys.path.insert(0,A+'/void_mass_network');import vmn_gen as G
web=G.generate(Mmin=3e15,seed=2000,comp_mode='lcdm');pp=web['params']
vv=np.linalg.norm(web['S'],axis=1)<pp['R_out'];nv=int(vv.sum());nc=len(web['MC'])
MSUN=1.989e30;Dm=-web['D'][vv]*MSUN;Wm=web['W'][vv]*MSUN;Cm=web['MC']*MSUN
CC=[-0.041,-0.035,0.1353,0.0];Rb=1.5*(1+1/2.937);aC=np.array([1,0,np.tan(np.radians(15))]);aC/=np.linalg.norm(aC);PC=42.0*aC
shp=lambda x:np.exp(sum(c*(np.minimum(np.abs(x),Rb)-1)**(k+1) for k,c in enumerate(CC)))
R0,R1=pp['gap_us']/1000,pp['R_out']/1000
xs=np.linspace((42-R1)/42,(42+R1)/42,2000);pmax=shp(xs).max()
def draw(n,rng):
  out=[]
  while sum(len(o) for o in out)<n:
    m=4*n;u=rng.random(m);r=(R0**3+u*(R1**3-R0**3))**(1/3);v=rng.normal(size=(m,3));v/=np.linalg.norm(v,axis=1)[:,None];P=r[:,None]*v
    out.append(P[rng.random(m)<shp(np.linalg.norm(P-PC,axis=1)/42)/pmax])
  return np.vstack(out)[:n]
for s in (1,2,3):
  rng=np.random.default_rng(s);PU=draw(nv+nc,rng)
  X=np.vstack([PU[:nv],PU[:nv],PU[nv:]]);M=np.concatenate([Dm,Wm,Cm])
  np.savez(f'densweb_s{s}.npz',X=X,M=M,nv=nv,nc=nc);print(s,len(M),M.sum(),flush=True)
np.savez('orig.npz',*G.kernel_bodies(web))
