import numpy as np,pickle,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from numpy.polynomial import legendre as Lg
o=np.load('rerun.npz')['o'];P=pickle.load(open('rerun.npz.paths.pkl','rb'))
psi=o[:,0];mu=np.cos(np.radians(psi));lz=np.log(o[:,1])
avg=lambda y:0.5*np.trapezoid(y[::-1],mu[::-1])
def fit(y):c=Lg.legfit(mu,y,6);return c,float(np.sqrt(np.mean((y-Lg.legval(mu,c))**2)))
S={};
for k in ('L','R','M'):
  v=[]
  for p in P:
    tot=p['Hr']*p['c2']+p['Hp']*p['s2'];f=(p['Hr']-p['Hp'])*(p['c2']-1/3)/tot;dl=np.diff(p['lnE'])
    fe={'L':f[:-1],'R':f[1:],'M':0.5*(f[:-1]+f[1:])}[k];v.append(np.sum(fe*dl))
  S[k]=np.array(v)
lnE_end=np.array([p['lnE'][-1] for p in P])
d=-(lz-avg(lz));st=-(S['M']-avg(S['M']));em=d-st
c,r=fit(d);cs,rs=fit(st);ce,re=fit(em)
out=dict(run='runs_sl2/slc_0.041037_0.035.npz, re-traced identically (dtt/rerun.npz; 1+z matches saved run exactly)',
 ratio=float(np.load('rerun.npz')['ratio']),lnE_path_vs_lnz_max=float(np.max(np.abs(lnE_end-lz))),
 a_l_total=c.tolist(),rms_resid_total=r,max_resid_total=float(np.max(np.abs(d-Lg.legval(mu,c)))),rms_DTT=float(np.sqrt(np.mean(d**2))),
 a_l_stretch=cs.tolist(),rms_resid_stretch=rs,a_l_emission_meanH=ce.tolist(),rms_resid_emission=re,
 stretch_sensitivity_left_vs_right_endpoint_max=float(np.max(np.abs((S['L']-avg(S['L']))-(S['R']-avg(S['R']))))),
 stretch_raw_range=[float(S['M'].min()),float(S['M'].max())],te_range=float(o[:,2].max()-o[:,2].min()),w_e_psi0_psi180=[float(o[0,4]),float(o[-1,4])])
json.dump(out,open('dtt_fit.json','w'),indent=1)
for k,v in out.items(): print(k,v)
ff=np.linspace(0,180,361);mf=np.cos(np.radians(ff));fig,ax=plt.subplots(figsize=(6,3.6))
ax.plot(180-psi,d*1e3,'ko',ms=3,label='ray-traced ΔT/T')
ax.plot(180-ff,Lg.legval(mf,c)*1e3,'k-',lw=1,label='Σ a_l P_l, l≤6')
ax.plot(180-psi,st*1e3,'r--',lw=1,label='stretch term');ax.plot(180-psi,em*1e3,'b:',lw=1.2,label='mean-H / emission-time term')
ax.set_xlabel('angle from CoE direction [deg]  (CoE left, ψ=0 right)');ax.set_ylabel('ΔT/T [10⁻³]');ax.legend(fontsize=7);ax.grid(alpha=.3)
plt.tight_layout();plt.savefig('dtt_fit.png',dpi=110)
