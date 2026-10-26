# Phase2 summary figure. Args: P5_NPZ(field run) P6_TAG. Panels: background, residual (bg+bodies), per-l rms, Migkas D profile, TT low-l, TT full.
import sys,json,glob,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase2/';p5c,t6=sys.argv[1],sys.argv[2]
A=np.load(O+'p5_R0.18957_s2000_L12_M1e+14_f1.npz');Ja=json.load(open(O+'p5_R0.18957_s2000_L12_M1e+14_f1.json'));C=np.load(O+p5c);Z=np.load(O+f'p6_{t6}.npz');J6=json.load(open(O+f'p6_{t6}.json'))
fig=plt.figure(figsize=(13,10))
hp.mollview(A['B']*1e6,sub=(3,2,1),title=f"CoE background D=34.96, rms {A['B'].std():.2e}",unit='ΔT/T ×1e6',fig=fig.number)
hp.mollview((A['Resid']-A['Resid'].mean())*1e6,sub=(3,2,2),title=f"+ bodies at R=0.190 Gpc: rms {Ja['rms_resid']:.1e}, max {Ja['max_resid']:.1e}",unit='ΔT/T ×1e6',fig=fig.number)
a=fig.add_subplot(3,2,3);L=np.arange(1,41);a.semilogy(L,Ja['amp_bg'][1:41],'o-',ms=3,label='background');a.semilogy(L,np.maximum(Ja['amp_resid'][1:41],1e-13),'s-',ms=3,label='bg + bodies');a.axhline(1e-5,color='r',ls='--',label='1e-5');a.set_xlabel('ℓ');a.set_ylabel('rms ΔT/T per ℓ');a.legend(fontsize=7)
a=fig.add_subplot(3,2,4);Ds,X2,R2=[],[],[]
for f in sorted(glob.glob(O+'p1_migkas_D*.json'),key=lambda f:float(f.split('_D')[1][:-5])):
    j=list(json.load(open(f)).values())[0];Ds.append(float(f.split('_D')[1][:-5]));X2.append(j['chi2']);R2.append(j['r2_pearson'])
a.plot(Ds,X2,'ko-');a.set_xlabel('D (Gpc)');a.set_ylabel('Migkas χ² (refit)');a.axvline(34.96,color='g',ls=':');b=a.twinx();b.plot(Ds,R2,'bs--');b.set_ylabel('r²',color='b');a.set_title('Migkas fit vs D (chosen D=34.96)')
Pf,Lt=Z['Pf'],Z['Lt'];ll=np.arange(len(Z['D']));x=J6['mean']['l2_2500']['chi2_planck'];y=J6['seed42']['l2_2500']['chi2_planck'];yc=J6['seed42_comp']['l2_2500']['chi2_planck']
for k,(lo,hi,xs) in enumerate([(2,30,'linear'),(2,2500,'log')]):
    a=fig.add_subplot(3,2,5+k);s=(Pf[:,0]>=lo)&(Pf[:,0]<=hi);a.errorbar(Pf[s,0],Pf[s,1],[Pf[s,2],Pf[s,3]],fmt='o',color='grey',ms=2,alpha=0.4,zorder=0,label='Planck');r=slice(lo,hi+1)
    a.plot(ll[r],Lt[r],'k',label='ΛCDM');a.plot(ll[r],Z['D'][r],'r',label=f'CoE ensemble χ²={x:.0f}');a.plot(ll[r],Z['Dr'][r],'c',lw=.7,label=f'seed 42 χ²={y:.0f}');a.plot(ll[r],Z['Drc'][r],'m--',lw=.8,label=f'seed 42 + bodies χ²={yc:.0f}');a.set_xscale(xs);a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.legend(fontsize=6)
fig.subplots_adjust(hspace=0.45,wspace=0.3,left=0.07,right=0.93,top=0.96,bottom=0.06);fig.savefig(O+'p7_summary.png',dpi=110);print('ok')
