import json,numpy as np,healpy as hp
exec(open('plot4.py').read())
z=np.load('jump_data.npz');w0,W,u,X,M,gid,sm=[z[k] for k in('w0','W','u','X','M','gid','sm')]
T0=2.7255e6;ll=np.arange(31);B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/reproduce_cmb'
def Dl(x):return l_(x)
def l_(x):return ll*(ll+1)*hp.anafast(np.asarray(hp.remove_dipole(x)),lmax=30)*T0**2/(2*np.pi)
P=np.loadtxt(B+'/planckTT.txt')
def chi(D):
  c=0
  for r in P:
    L=int(r[0])
    if 2<=L<=30:e=r[3] if D[L]>r[1] else r[2];c+=((D[L]-r[1])/e)**2
  return c
base=sm+w0;print('check 1.5%',chi(l_(base)),' 1.75%',chi(l_(base+W.sum(0))))
F=np.round(np.arange(0.015,0.01750001,0.0001),5);cs=np.cumsum(W,0);rows=[]
for f in F:
  n=int((u<f).sum());m=base+(cs[n-1] if n else 0);D=l_(m);rows.append(dict(f=float(f),n_added=n,chi2=chi(D),r2=r2(D),Dl=D.tolist()))
for r in rows:print('%.4f n=%3d chi2=%.1f r2=%.3f D20=%.0f D25=%.0f D30=%.0f'%(r['f'],r['n_added'],r['chi2'],r['r2'],r['Dl'][20],r['Dl'][25],r['Dl'][30]))
# per-body: chi2 change adding body alone to 1.5% map, and removing from 1.75% map
full=base+W.sum(0);c0=chi(l_(base));c1=chi(l_(full))
hi=lambda D:np.mean(D[20:31])
per=[]
for i in range(len(u)):
  Da=l_(base+W[i]);Dr=l_(full-W[i])
  per.append(dict(i=i,u=float(u[i]),type='void' if gid[i]<41000000 else '',M=float(M[i]),X=X[i].tolist(),r=float(np.linalg.norm(X[i])),dchi_add=chi(Da)-c0,dchi_remove=chi(Dr)-c1,dhi_add=hi(Da)-hi(l_(base)),dhi_rem=hi(Dr)-hi(l_(full))))
json.dump(dict(rows=rows,per=per),open('thin4_jump.json','w'))
o=sorted(per,key=lambda p:p['dhi_rem'])[:8]
for p in o:print(p['i'],'%.5f M=%.2e r=%.3f |r-4|=%.4f dchi_add=%.1f dchi_rem=%.1f dhi_add=%.0f dhi_rem=%.0f'%(p['u'],p['M'],p['r'],abs(p['r']-4),p['dchi_add'],p['dchi_remove'],p['dhi_add'],p['dhi_rem']))
