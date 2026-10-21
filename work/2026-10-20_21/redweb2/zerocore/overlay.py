import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt,zc
d=json.load(open('scanB.json'))+json.load(open('scan.json'));ok=[o for o in d if 'chi2' in o]
RM=1.92  # deepest Migkas cluster, z=0.447 -> cz/H0
A=min([o for o in ok if o['D']>RM],key=lambda o:o['chi2']);Bu=min(ok,key=lambda o:o['chi2']);Bi=min([o for o in ok if o['D']<=RM],key=lambda o:o['chi2'])
red=json.load(open('/workspace/redweb2/densweb/migkas_densweb_s1.json'))['smooth']['plane_fit']
o=np.argsort(zc.L);fig,ax=plt.subplots(2,1,figsize=(12,13))
a=ax[0];a.errorbar(zc.L[o],zc.Y[o],yerr=zc.Sg[o],fmt='o',c='k',ms=4,capsize=2,label='Migkas 2020, 71 bins')
a.axhline(1,c='0.5',ls=':',label=f'flat sky: χ²={zc.chi(np.ones(71)):.1f}')
a.plot(zc.L[o],np.array(red['A'])[o],'r-',lw=2,label=f"red (poly ln ρ, D=42): χ²={red['chi2']:.1f}, Pearson r²={red['pearson_r2']:.2f}")
cs={'A':'tab:blue','B':'tab:green'};runs=[('A (D>R=1.92 Gpc)',A,'tab:blue','-')]
runs.append(('B unrestricted',Bu,'tab:green','--') if Bu is not A else ('B best with D≤R',Bi,'tab:green','--'))
for lab,r,c,ls in runs:
    a.plot(zc.L[o],np.array(r['A'])[o],c=c,ls=ls,lw=2,label=f"zero-CoE {lab}: D={r['D']:g} Gpc, a={r['a']}, b={r['b']}, r_pk={r['xm']}D, frac={r['frac']:g}: χ²={r['chi2']:.1f}, Pearson r²={r['pearson_r2']:.2f}")
a.set_xlabel('l (deg)');a.set_ylabel('A/A_all');a.legend(loc='lower center',fontsize=8.5);a.set_title('Migkas A/A_all, zero-CoE-density law ρ ∝ (r/r_p)^a exp(−(r/r_p)^b); real z, no vp, no Λ')
b=ax[1];RR=14.3
for lab,r,c,ls in runs:
    D=r['D'];rr=np.linspace(0,D+RR,800);xp=r['xm']/(r['a']/r['b'])**(1/r['b']);y=zc.shape_f(r['a'],r['b'],xp)(rr/D)*((rr/D)<=1.5*(1+RR/D))
    b.plot(rr,y/y.max(),c=c,ls=ls,lw=2,label=f'zero-CoE {lab} (D={D:g})');b.axvline(D,c=c,ls=':')
rr=np.linspace(0,42+RR,800);CC=[-0.041,-0.035,0.1353,0];Rb=1.5*(1+1/2.937);x=np.minimum(rr/42,Rb)
yr=np.exp(sum(c*(x-1)**(k+1) for k,c in enumerate(CC)));b.plot(rr,yr/yr.max(),'r',lw=2,label='red ρ (D=42)')
for v in (42-RR,42,42+RR):b.axvline(v,c='r',ls=':',lw=.8)
b.set_xlabel('r from CoE (Gpc)');b.set_ylabel('ρ / max');b.legend(loc='lower center',fontsize=9);b.set_xlim(0,42+RR);b.set_title('ρ(r); dotted = D of each model, red dotted = 42∓14.3')
plt.tight_layout();plt.savefig('migkas_zero_overlay.png',dpi=100)
J=dict(R_cut_Gpc=RM,red=dict(chi2=red['chi2'],pearson_r2=red['pearson_r2'],CoE_l=red['CoE_l']),flat_chi2=zc.chi(np.ones(71)),
 A={k:A[k] for k in A if k!='A'},B_unrestricted={k:Bu[k] for k in Bu if k!='A'},B_best_D_le_R={k:Bi[k] for k in Bi if k!='A'},n_ok=len(ok),n_tried=len(d))
json.dump(J,open('migkas_zero_fits.json','w'),indent=1)
open('READY_migkas','w').write('/workspace/redweb2/zerocore/migkas_zero_overlay.png\n');print(json.dumps(J,indent=0)[:1500])
