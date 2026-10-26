# rev6 Fig 18 (P_delta split into terms; radial window W(chi)) regenerated from the rev7 field fit p8_field_consistent.json
# (D=34.96 Gpc, R=0.18957 Gpc). Same functional form as Pd() in joint_solve/phase2/p8_field_consistent.py; W = Gaussian(chiw, sw) as in acoustic_model.py.
import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
J=json.load(open('/workspace/coepdf/joint_solve/phase2/p8_field_consistent.json'));p=J['best_p'];R=J['R_Gpc'];s=R/4;K0=100/s
k=np.logspace(np.log10(1e-2/s),np.log10(200/s),2000);lk=np.log(k/K0)
base=10**p['lgA']*k**(p['n']+p['alpha']*lk);damp=np.exp(-(k/p['kD'])**p['m'])
ac=p['Aw']*np.exp(p['beta']*lk)*np.cos(k*p['rs']+p['phi']);hs=p['B']*np.exp(-k/p['kB'])*np.cos(k*p['r2']+p['phi2'])
fig,(a,b)=plt.subplots(1,2,figsize=(10,3.8))
a.loglog(k,base*damp*(1+ac+hs),'k',label='total P_δ');a.loglog(k,base*damp,'b--',label='smooth × damping')
a.loglog(k,np.abs(base*damp*ac),'r',lw=.7,label='|acoustic term|');a.loglog(k,np.abs(base*damp*hs),'g',lw=.7,label='|half-scale term|')
a.set_xlabel('k (Gpc⁻¹)');a.set_ylabel('P_δ (arb.)');a.set_ylim(base[-1]*damp[-1]*1e-3+1e-300,None);a.legend(fontsize=7)
ch=np.linspace(0,R,1000);b.plot(ch,np.exp(-0.5*((ch-p['chiw'])/p['sw'])**2),'k');b.axvline(R,color='r',ls=':',label=f'R = {R} Gpc')
b.set_xlabel('χ (Gpc)');b.set_ylabel('W(χ)');b.set_title(f"χ_w={p['chiw']:.4f}, s_w={p['sw']:.4f} Gpc",fontsize=9);b.legend(fontsize=7)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=130);print('ok')
