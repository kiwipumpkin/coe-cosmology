# rev7 regeneration of rev6 Fig 6 (Migkas direct fits, ray-traced CMB, density shape, low multipoles). rev6 (poly_plots.py) drew the
# D=42 Gpc runs. This uses the rev7 solve products: p1_migkas_D*.json (profile refits vs D), p2_trace_D*.npz (famP.py ray traces,
# same shape law), p2_sm_D*_ns64.npy (background skies), with D=34.96 Gpc the adopted value.
import json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/coepdf/joint_solve/phase2/';d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True)
l,y=d['l_deg'],d['A_over_Aall'];el,eh=d['A_err_low'],d['A_err_high']
PL=np.load(O+'p6_p8.npz')['Pf'];CC=[-0.041,-0.035,0.1353,0];RV=14.30048287260589
Ds=['30','34.96','42','50','58.2'];col=dict(zip(Ds,['C1','C3','C0','C2','C4']));out={}
fig,ax=plt.subplots(1,4,figsize=(25,5.6))
a=ax[0];a.errorbar(l,y,[el,eh],fmt='o',ms=3,color='k',alpha=.6,label='Migkas 71 cones')
for D in Ds:
    p=list(json.load(open(O+f'p1_migkas_D{D}.json')).values())[0];c=np.array(p['curve']);lf=np.arange(len(c))*360/len(c);out[D]=dict(chi2=p['chi2'],r2=p['r2_pearson'],l0=p['l0'])
    a.plot(lf,c,color=col[D],lw=2.4 if D=='34.96' else 1,ls='-' if D=='34.96' else '--',label=f'D={D} Gpc: χ²={p["chi2"]:.2f}, r²={p["r2_pearson"]:.3f}, ℓ0={p["l0"]:.1f}°')
a.set_xlim(0,360);a.set_xlabel('galactic longitude [deg]');a.set_ylabel('A/A_all');a.legend(fontsize=7);a.set_title('Migkas A/A_all, profile refit vs D')
for D in ('34.96','42','58.2'):
    o=np.load(O+f'p2_trace_D{D}.npz')['o'];o=o[np.argsort(o[:,0])];psi=o[:,0];lz=np.log(o[:,1]);w=np.sin(np.radians(psi))
    dT=-(lz-np.trapezoid(lz*w,psi)/np.trapezoid(w,psi));ax[1].plot(psi,dT*1e5,color=col[D],marker='o',ms=2,lw=2 if D=='34.96' else 1,label=f'D={D}')
    sm=np.load(O+f'p2_sm_D{D}_ns64.npy');cl=hp.anafast(sm,lmax=6);T0=2.7255e6;Dl=np.arange(7)*(np.arange(7)+1)*cl/2/np.pi*T0**2;out[D]['Dl2_6']=Dl[2:7].tolist();out[D]['rms']=float(sm.std())
    x=np.linspace(1e-3,2.6,2000);Rb=1.5*(1+RV/float(D));ax[2].plot(x,np.exp(sum(c*(np.minimum(x,Rb)-1)**(k+1) for k,c in enumerate(CC)))*(x<=Rb),color=col[D],label=f'D={D}, R_b={Rb:.2f}')
ax[1].invert_xaxis();ax[1].set_xlabel('ψ [deg] (180 = toward CoE)');ax[1].set_ylabel('ΔT/T [10⁻⁵] (monopole removed)');ax[1].set_title('Ray-traced background CMB (famP poly)');ax[1].legend(fontsize=7)
ax[2].axvline(1,color='k',ls=':');ax[2].set_xlabel('distance from CoE / D (observer=1)');ax[2].set_ylabel('initial ρ shape');ax[2].set_title('Density shape (same c_k; R_b depends on D)');ax[2].legend(fontsize=7)
a=ax[3];k=np.arange(2,7)
for i,D in enumerate(('34.96','42','58.2')):a.bar(k-0.27+0.27*i,out[D]['Dl2_6'],0.27,color=col[D],label=f'background sky D={D}')
s=PL[:,0]<=6;a.errorbar(PL[s,0],PL[s,1],[PL[s,2],PL[s,3]],fmt='k_',ms=20,mew=2,capsize=4,label='Planck 2018 TT')
a.set_yscale('log');a.set_xlabel('ℓ');a.set_ylabel('D_ℓ [μK²] (before bodies)');a.set_title('Background multipoles vs Planck (bodies remove these, §CMB)');a.legend(fontsize=7)
fig.tight_layout();fig.savefig(__file__.replace('.py','.png'),dpi=90);json.dump(out,open(__file__.replace('.py','.json'),'w'),indent=1);print(out)
