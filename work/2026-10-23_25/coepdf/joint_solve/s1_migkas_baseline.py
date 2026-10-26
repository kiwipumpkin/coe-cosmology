# Step 0: Migkas A/A_all baseline (locked). Recompute chi2 and Pearson r2 of existing CatWISE-locked best curves vs 71 digitised bins.
import numpy as np,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
N='/workspace/cosmology-course/notes/'
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
l,y=d['l_deg'],d['A_over_Aall'];s=0.5*(d['A_err_low']+d['A_err_high'])
cur={'eta_edge_free (D=34.96, l0=19.9)':N+'catwise-eta-migkas-fit-d14.3.csv','eta0.25 l0=15.2':N+'catwise-density-levelset-migkas-eta_fixed15.2.csv',
     'soft l0=15.2':N+'catwise-density-levelset-migkas-soft_fixed15.2.csv','exp l0=15.2':N+'catwise-density-levelset-migkas-exp_fixed15.2.csv'}
out={'chi2_null':float(np.sum((y-1)**2/s**2)),'errors':'mean of asymmetric Migkas errors','fits':{}}
fig,a=plt.subplots(figsize=(9,4.5));a.errorbar(l,y,s,fmt='ko',ms=3,label='Migkas A/A_all (71 bins)')
for k,f in cur.items():
    c=np.loadtxt(f,delimiter=',',skiprows=1);m=np.interp(l,c[:,0],c[:,1],period=360)
    out['fits'][k]=dict(chi2=float(np.sum((y-m)**2/s**2)),r2=float(np.corrcoef(y,m)[0,1]**2),src=f);a.plot(c[:,0],c[:,1],label=f"{k}: χ²={out['fits'][k]['chi2']:.2f}, r²={out['fits'][k]['r2']:.3f}")
a.set_xlabel('Galactic longitude (deg)');a.set_ylabel('A/A_all');a.legend(fontsize=7);fig.tight_layout();fig.savefig('s1_migkas_baseline.png',dpi=110)
json.dump(out,open('s1_migkas_baseline.json','w'),indent=1);print(json.dumps(out,indent=1))
