import glob,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from numpy.polynomial.legendre import legval
from scipy.interpolate import CubicSpline,interp1d
x=np.linspace(-1,1,8001);ps=np.degrees(np.arccos(x));PL=0.356;PD=71.2  # Planck per-l rms (1e-5): l=2, l=1 (3362uK/sqrt3)
def alm(f):
    o=np.load(f)['o'];o=o[np.argsort(o[:,0])];y=-CubicSpline(o[:,0],np.log(o[:,1]))(ps);y-=y.mean();r=[]
    for l in range(1,6):c=np.zeros(l+1);c[l]=1;r.append((2*l+1)/2*np.trapezoid(y*legval(x,c),x)/np.sqrt(2*l+1)/1e-5)
    return r,float(np.mean(o[:,1]))
MP=json.load(open('migkas_vs_ratio.json'));rr=np.array(sorted(map(float,MP)));cc=np.array([MP[str(k) if str(k) in MP else repr(k)]['chi2'] for k in rr]);R2=np.array([MP[str(k)]['R2'] for k in rr])
chif=interp1d(rr,cc,bounds_error=False);r2f=interp1d(rr,R2,bounds_error=False);cmin=cc.min()
rows=[];fams={}
for f in sorted(glob.glob('runs/*.npz')):
    fam,p,m=f[5:-4].split('_');m=float(m);key=f'{fam} {p}'
    try:a,z=alm(f);d=np.load(f);rows.append(dict(family=key,frac=m,ratio=float(d['ratio']),zmean=z,a=a))
    except Exception as e:print('fail',f,e)
out=[]
for r in rows:
    k0=[q for q in rows if q['family']==r['family'] and q['frac']==0]
    kin=k0[0]['a'] if k0 else [np.nan]*5;grav=[t-k for t,k in zip(r['a'],kin)]
    c=float(chif(r['ratio']));e=dict(family=r['family'],M_over_Mcrit=r['frac'],ratio=r['ratio'],mean_1pz=r['zmean'],migkas_chi2=c,migkas_dchi2=c-cmin,migkas_R2=float(r2f(r['ratio'])),
      net_l1=r['a'][0],net_l2=r['a'][1],net_l3=r['a'][2],net_l4=r['a'][3],grav_l2=grav[1],kin_l2=kin[1],
      cancel_frac=(-kin[1]/grav[1] if r['frac']>0 and grav[1] else np.nan),l2_over_planck=abs(r['a'][1])/PL,dipole_xPlanck=abs(r['a'][0])/PD,v_needed_kms=370*abs(r['a'][0])/PD)
    out.append(e)
json.dump(out,open('scan_results.json','w'),indent=1)
print(f"{'family':9s} {'M/Mc':>5} {'ratio':>6} {'1+z':>6} {'dchi2':>6} {'R2':>5} {'l1':>8} {'l2':>8} {'grav2':>8} {'kin2':>7} {'canc':>5} {'l3':>7} {'l4':>6} {'v_kms':>7}")
for e in sorted(out,key=lambda e:(e['family'],e['M_over_Mcrit'])):
    print(f"{e['family']:9s} {e['M_over_Mcrit']:5.2f} {e['ratio']:6.3f} {e['mean_1pz']:6.0f} {e['migkas_dchi2']:6.2f} {e['migkas_R2']:5.3f} {e['net_l1']:8.1f} {e['net_l2']:8.1f} {e['grav_l2']:8.1f} {e['kin_l2']:7.1f} {e['cancel_frac']:5.2f} {e['net_l3']:7.1f} {e['net_l4']:6.1f} {e['v_needed_kms']:7.0f}")
fig,ax=plt.subplots(1,2,figsize=(17,7))
for fam in sorted(set(e['family'] for e in out)):
    E=sorted([e for e in out if e['family']==fam],key=lambda e:e['M_over_Mcrit'])
    ax[0].plot([e['M_over_Mcrit'] for e in E],[abs(e['net_l2']) for e in E],'o-',label=fam)
    ax[1].plot([e['ratio'] for e in E],[abs(e['net_l2']) for e in E],'o-',label=fam)
ok=rr[cc-cmin<4];ax[1].axvspan(ok.min(),ok.max(),color='g',alpha=.12,label='Migkas Δχ²<4');ok1=rr[cc-cmin<1];ax[1].axvspan(ok1.min(),ok1.max(),color='g',alpha=.2,label='Migkas Δχ²<1')
for a in ax:a.axhline(PL,c='gold',lw=3,label='Planck l=2 (9.7 μK)');a.set_yscale('log');a.set_ylabel('|net l=2| (1e-5)');a.legend(fontsize=9)
ax[0].set_xlabel('mass M/M_crit inside our radius (4.36 = Migkas-fitted "b")');ax[1].set_xlabel('H_par/H_perp at us (sets Migkas fit)')
fig.suptitle('Net CMB quadrupole vs mass, per density family (all re-locked to CatWISE log-slope 0.041), 1+z≈1100, shrinking radius')
plt.savefig('l2_vs_mass.png',dpi=110,bbox_inches='tight')
