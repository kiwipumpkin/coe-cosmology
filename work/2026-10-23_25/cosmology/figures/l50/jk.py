import numpy as np,sys,csv,json
sys.path.insert(0,'/workspace/cosmology/figures')
from make_A_vs_angle_axisymmetric_before_rotation import periodic_interp
L=[q for q in open('/workspace/cosmology/migkas/migkas2020_tablec1_vizier.tsv') if not q.startswith('#') and q.strip()]
h=L[0].rstrip('\n').split('\t');R=[r.rstrip('\n').split('\t') for r in L[3:]]
g_=lambda k,f=float:np.array([f(r[h.index(k)]) if f is str or r[h.index(k)].strip() else np.nan for r in R])
name=g_('Cluster',str);z=g_('z');l=g_('GLON');b=g_('GLAT');T=g_('T');LX=g_('LX');NH=g_('NHtot')
E=np.sqrt(.3*(1+z)**3+.7);y=np.log10(LX/100/E);x=np.log10(T/4)
grid=np.arange(5,360,5.)
def fit(m):
  Bs=np.cov(x[m],y[m],bias=True)[0,1]/np.var(x[m]);return y[m].mean()-Bs*x[m].mean(),Bs
def curves(keep):
  la,Ball=fit(keep);A=[];B=[]
  for g in grid:
    m=keep&(np.abs((l-g+180)%360-180)<=45);a,bb=fit(m);A.append(10**(a-la));B.append(bb/Ball)
  return np.array(A),np.array(B)
dm=np.genfromtxt('/workspace/cosmology/map_compare/amp_mean_param_noise/A_vs_angle_ideal_plus_chris_rotation.csv',delimiter=',',names=True)
AFm=periodic_interp(dm['ell_deg'],dm['AF_smooth'],grid)
gb=np.genfromtxt('/workspace/cosmology/figures/fourier_resid/sheet_BL_copy2_density_w_rotat.csv',delimiter=',',names=True)
BLm=periodic_interp(gb['ell_deg'],gb['BL_modelB_sheet'],grid)
allk=np.ones(len(z),bool);A0,B0=curves(allk)
la,Bg=fit(allk);res=y-(la+Bg*x)  # dex residual from global relation
S=[np.where(grid==c)[0][0] for c in (40,50,60)];Sw=[np.where(grid==c)[0][0] for c in range(30,75,5)]
print('full rA(40,50,60)',(A0-AFm)[S].round(3),'rB',(B0-BLm)[S].round(3))
inw=np.abs((l-50+180)%360-180)<=65  # in any slice centred 30..70
idx=np.where(inw)[0];rows=[]
for i in idx:
  k=allk.copy();k[i]=False;A1,B1=curves(k)
  dA=(A1-A0)[S];dB=(B1-B0)[S]
  rows.append(dict(i=int(i),name=name[i].strip(),l=l[i],b=b[i],z=z[i],T=T[i],LX=LX[i],NH=NH[i],res=res[i],dA=dA.mean(),dB=dB.mean(),
   bumpA=((A1-AFm)[S]).mean(),bumpB=((B1-BLm)[S]).mean()))
bA0=(A0-AFm)[S].mean();bB0=(B0-BLm)[S].mean()
for r in rows:r['infl']=abs(r['dA'])/max(abs(bA0),1e-9)+abs(r['dB'])/max(abs(bB0),1e-9)
rows.sort(key=lambda r:-(-r['dA']*np.sign(bA0)-r['dB']*np.sign(bB0)*abs(bA0/bB0)))
print('N in window',len(idx),'bumpA',round(bA0,4),'bumpB',round(bB0,4))
json.dump(dict(rows=rows,bA0=bA0,bB0=bB0),open('jk.json','w'),default=float)
w=csv.writer(open('clusters_l30_70.csv','w'));w.writerow(['name','l','b','z','T_keV','LX_1e44','NHtot_1e20','LXT_resid_dex','dA_jk','dB_jk','bumpA_after','bumpB_after'])
for r in sorted(rows,key=lambda r:r['l']):w.writerow([r['name'],r['l'],r['b'],r['z'],r['T'],r['LX']/100,r['NH'],round(r['res'],3),round(r['dA'],4),round(r['dB'],4),round(r['bumpA'],4),round(r['bumpB'],4)])
print('rank name l b z T LX44 NH res dA dB')
for k,r in enumerate(rows[:15]):print(k+1,r['name'],f"{r['l']:.1f} {r['b']:+.1f} {r['z']:.3f} {r['T']:.1f} {r['LX']/100:.2f} {r['NH']:.1f} {r['res']:+.2f} {r['dA']:+.4f} {r['dB']:+.4f}")
# top-k removal
top=[r['i'] for r in rows]
out={'grid':grid.tolist(),'A0':A0.tolist(),'B0':B0.tolist(),'AF':AFm.tolist(),'BL':BLm.tolist()}
for kk in (3,5,10):
  k=allk.copy();k[top[:kk]]=False;A1,B1=curves(k);out[f'A{kk}']=A1.tolist();out[f'B{kk}']=B1.tolist()
  print('top',kk,'bumpA',round((A1-AFm)[S].mean(),4),'bumpB',round((B1-BLm)[S].mean(),4),'chiwin-ish rmsA30-70',round(np.sqrt(np.mean((A1-AFm)[Sw]**2)),4))
# random control: remove kk random window clusters
rng=np.random.default_rng(0)
for kk in (5,10):
  v=[]
  for _ in range(200):
    k=allk.copy();k[rng.choice(idx,kk,replace=False)]=False;A1,B1=curves(k);v.append(((A1-AFm)[S].mean(),(B1-BLm)[S].mean()))
  v=np.array(v);print('random',kk,'bumpA mean/sd',v[:,0].mean().round(4),v[:,0].std().round(4),'bumpB',v[:,1].mean().round(4),v[:,1].std().round(4))
json.dump(out,open('curves.json','w'))
