import numpy as np,json,os,sys,csv
os.chdir('/workspace/cosmology/figures/fourier_resid');sys.path.insert(0,'/workspace/cosmology/figures');sys.path.insert(0,'.')
exec(open('recon_n3.py').read().split("cols=['purple'")[0])
lg=l.copy();Th=T
J=open('../l50/jk.py').read().split("allk=np.ones")[0]
exec(J)  # defines cluster arrays x,y,l (overwrites l!), grid, fit, AFm,BLm
lc=l.copy();l=lg
def chain(Ag,Bg_):
  A_=periodic_interp(grid,Ag,l);B_=periodic_interp(grid,Bg_,l)
  rA=A_-AFl;rB=B_-BLl;st=[]
  for n in range(1,4):
    h=np.column_stack([np.cos(n*th),np.sin(n*th)]);X=np.vstack([np.column_stack([o,Z,h]),np.column_stack([Z,o,h])]);s=np.r_[sA,sB]
    p=np.linalg.lstsq(X/s[:,None],np.r_[rA,rB]/s,rcond=None)[0];st.append((n,p));rA=rA-X[:71]@p;rB=rB-X[71:]@p
  out={}
  for nm,Y,mb,k in (('A',A_,AFl,0),('B',B_,BLl,1)):
    m=mb.copy()
    for n,p in st:m=m+Th(n,p,th,k)
    r=100*(Y-m)/Y
    out[nm]={w:np.array([np.mean(r[np.abs((l-c+180)%360-180)<=w]) for c in l]) for w in (15,25,35)}
  return out
def curves_idx(ix):
  xx,yy,ll=x[ix],y[ix],lc[ix]
  f=lambda m:(yy[m].mean()-(np.cov(xx[m],yy[m],bias=True)[0,1]/np.var(xx[m]))*xx[m].mean(),np.cov(xx[m],yy[m],bias=True)[0,1]/np.var(xx[m]))
  a0,b0=f(np.ones(len(ix),bool));A=[];B=[]
  for gg in grid:
    m=np.abs((ll-gg+180)%360-180)<=45;a,bb=f(m);A.append(10**(a-a0));B.append(bb/b0)
  return np.array(A),np.array(B)
N=len(x);base=chain(*curves_idx(np.arange(N)))
M=json.load(open('/tmp/M25.json'))
dA=np.max(np.abs(base['A'][25]-np.array(M['A'])));dB=np.max(np.abs(base['B'][25]-np.array(M['B'])))
print('N clusters',N,'repro max dev A',dA,'B',dB,flush=True)
spots=[('A','CMB dipole cold',97.2,-1),('A','Bulk flow away',113.9,-1),('A','CMB Cold Spot',156.2,-1),('A','CMB dipole hot',277.2,1),('A','Bulk flow toward',293.9,1),
('B','NVSS cold',58.5,-1),('B','Hemi asym cold',29.9,-1),('B','CatWISE cold',72.2,-1),('B','NVSS hot',238.5,1),('B','Hemi asym hot',209.9,1),('B','CatWISE hot',252.2,1),('B','Loop I',335.1,1)]
ev=lambda o,K,w,lm:float(periodic_interp(l,o[K][w],np.array([lm]))[0])
rng=np.random.default_rng(1);ND=2000;D=[]
for i in range(ND):
  D.append(chain(*curves_idx(rng.integers(0,N,N))))
rows=[]
for K,nm,lm,sg in spots:
  r={'K':K,'spot':nm,'l_m':lm,'need':'>0' if sg>0 else '<0'}
  for w in (15,25,35):
    v=ev(base,K,w,lm);bs=np.array([ev(d,K,w,lm) for d in D]);s=bs.std(ddof=1)
    r[f'val{w}']=round(v,3);r[f'sig{w}']=round(s,3);r[f'z{w}']=round(v/s,2);r[f'signfrac{w}']=round(np.mean(np.sign(bs)==sg),3)
  rows.append(r);print(r)
w=csv.DictWriter(open('/workspace/cosmology/figures/l50/spot_residual_errors.csv','w'),fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
