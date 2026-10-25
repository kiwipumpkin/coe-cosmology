import numpy as np,csv,sys,json,warnings;warnings.filterwarnings('ignore')
sys.path.insert(0,'/workspace/cosmology/figures')
from make_A_vs_angle_axisymmetric_before_rotation import periodic_interp
from scipy.integrate import quad
rng=np.random.default_rng(1)
z,l,b,T,LX=np.load('clusters.npy').T
E=lambda z:np.sqrt(.3*(1+z)**3+.7)
dc=lambda z:quad(lambda x:1/E(x),0,z)[0]*299792.458/70
y=np.log10(LX/100/E(z));x=np.log10(T/4)
grid=np.arange(5,360,5.)  # 71 centers matching Migkas bins
def curves(idx):
  Y,X,Lg=y[idx],x[idx],l[idx]
  def fit(m):
    Xm,Ym=X[m],Y[m];Bs=np.cov(Xm,Ym,bias=True)[0,1]/np.var(Xm);a=Ym.mean()-Bs*Xm.mean();return 10**a,Bs
  Aall,Ball=fit(np.ones(len(Y),bool));A=[];B=[];N=[]
  for g in grid:
    m=np.abs((Lg-g+180)%360-180)<=45;a,bb=fit(m);A.append(a/Aall);B.append(bb/Ball);N.append(m.sum())
  return np.array(A),np.array(B),np.array(N)
# models
dm=np.genfromtxt('/workspace/cosmology/map_compare/amp_mean_param_noise/A_vs_angle_ideal_plus_chris_rotation.csv',delimiter=',',names=True)
AFm=periodic_interp(dm['ell_deg'],dm['AF_smooth'],grid)
gb=np.genfromtxt('/workspace/cosmology/figures/fourier_resid/sheet_BL_copy2_density_w_rotat.csv',delimiter=',',names=True)
BLm=periodic_interp(gb['ell_deg'],gb['BL_modelB_sheet'],grid)
th=np.radians(grid)
def harm(r):
  X=np.column_stack([np.ones_like(th)]+[f(k*th) for k in range(1,5) for f in(np.cos,np.sin)])
  p=np.linalg.lstsq(X,r,rcond=None)[0];a,s=p[1::2],p[2::2];n=np.arange(1,5)
  return np.hypot(a,s),(np.degrees(np.arctan2(s,a))%360)/n
def circ(ph,n):  # ph in deg, period 360/n
  ang=np.radians(ph*n);m=np.angle(np.mean(np.exp(1j*ang)));R=np.abs(np.mean(np.exp(1j*ang)))
  return (np.degrees(m)%360)/n,np.degrees(np.sqrt(-2*np.log(max(R,1e-9))))/n
# validate full sample vs digitized
rows=list(csv.DictReader(open('/workspace/cosmology/migkas/fig6_digitized.csv')))
Ad=np.array([float(r['A_over_Aall']) for r in rows]);Bd=np.array([float(r['B_over_Ball']) for r in rows])
Af,Bf,Nf=curves(np.arange(len(z)))
print('validation full sample: r(A,digitized)=%.3f r(B,digitized)=%.3f; A range %.3f-%.3f'%(np.corrcoef(Af,Ad)[0,1],np.corrcoef(Bf,Bd)[0,1],Af.min(),Af.max()))
order=np.argsort(z);shells={'full':order}
for k,part in enumerate(np.array_split(order,3)):shells[f'z{k+1}']=part
out={}
for name,idx in shells.items():
  A,B,N=curves(idx);rA,rB=A-AFm,B-BLm
  CA,PA=harm(rA);CB,PB=harm(rB)
  bsA=[];bsB=[]
  for _ in range(500):
    s=rng.choice(idx,len(idx));a,bb,_n=curves(s)
    if not(np.all(np.isfinite(a)) and np.all(np.isfinite(bb))):continue
    bsA.append(harm(a-AFm));bsB.append(harm(bb-BLm))
  res={'N':len(idx),'zmed':float(np.median(z[idx])),'zmin':float(z[idx].min()),'zmax':float(z[idx].max())}
  res['d_Mpc']=dc(res['zmed']);res['minN_per_slice']=int(N.min())
  res['n_boot_ok']=len(bsA)
  for lab,C,P,bs in(('A',CA,PA,bsA),('B',CB,PB,bsB)):
    Cb=np.array([q[0] for q in bs]);Pb=np.array([q[1] for q in bs]);res[lab]=[]
    for i in range(4):
      n=i+1;mph,sph=circ(Pb[:,i],n)
      res[lab].append(dict(n=n,C=float(C[i]),C_err=float(Cb[:,i].std()),C_bias_med=float(np.median(Cb[:,i])),peak=float(P[i]),peak_err=float(sph),lam_Mpc=2*np.pi*res['d_Mpc']/n))
  res['curveA']=A.tolist();res['curveB']=B.tolist();out[name]=res
  print(name,res['N'],round(res['zmed'],4),round(res['d_Mpc']),'minN',res['minN_per_slice'],'boot',len(bsA))
  for lab in 'AB':
    for q in res[lab]:print('  ',lab,q['n'],'C=%.3f±%.3f peak=%.0f±%.0f lam=%.0f Mpc'%(q['C'],q['C_err'],q['peak'],q['peak_err'],q['lam_Mpc']))
json.dump(out,open('zshell_results.json','w'),indent=1)
