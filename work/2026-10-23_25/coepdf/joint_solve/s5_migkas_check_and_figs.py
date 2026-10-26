# Step 4: (i) Migkas check: perturbation of A/A_all by the compensating bodies, delta A/A ~ 2*dPhi/c^2 averaged along each
# Migkas cone direction (b=0, 90-deg cone approximated by sightline average over |dl|<=45, 0.5-14.3 Gpc); recompute chi2/r2 of the
# locked eta-edge curve. (ii) Figures: sky before/after, low-l dT/T, body positions, TT spectrum.
import sys,json,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
s2,s3,s4=sys.argv[1:4];Z=np.load(s2);J2=json.load(open(s2.replace('.npz','.json')));Z3=np.load(s3);J4=json.load(open(s4))
G=4.30e-12;c2=(2.998e5)**2;y,M=Z['y'],Z['M']
d=np.genfromtxt('/workspace/cosmology/migkas/fig6_digitized.csv',delimiter=',',names=True);l,yA=d['l_deg'],d['A_over_Aall'];s=0.5*(d['A_err_low']+d['A_err_high'])
c=np.loadtxt('/workspace/cosmology-course/notes/catwise-eta-migkas-fit-d14.3.csv',delimiter=',',skiprows=1);m0=np.interp(l,c[:,0],c[:,1],period=360)
rr=np.linspace(0.5,14.3,40);dPhi=[]
for L0 in l:
    v=[]
    for dl in np.linspace(-45,45,19):
        t=np.radians(L0+dl);u=np.array([np.cos(t),np.sin(t),0.]);P=rr[:,None]*u
        v.append(np.mean(-G*np.sum(M[None]/np.sqrt(((P[:,None,:]-y[None])**2).sum(-1)+0.01),1)/c2))
    dPhi.append(np.mean(v))
dPhi=np.array(dPhi);m1=m0*(1+2*(dPhi-dPhi.mean()))
mig=dict(dPhi_over_c2_range=[float(dPhi.min()),float(dPhi.max())],chi2_before=float(np.sum((yA-m0)**2/s**2)),chi2_after=float(np.sum((yA-m1)**2/s**2)),
 r2_before=float(np.corrcoef(yA,m0)[0,1]**2),r2_after=float(np.corrcoef(yA,m1)[0,1]**2),max_dA=float(np.abs(m1-m0).max()))
json.dump(mig,open('s5_migkas_check.json','w'),indent=1);print(mig)
fig=plt.figure(figsize=(12,9))
hp.mollview(Z['B']*1e6,sub=(3,2,1),title=f"CoE background (sm+rs_cc), rms {Z['B'].std():.2e}",unit='ΔT/T ×1e6',fig=fig.number)
hp.mollview((Z['Tot']-Z['Tot'].mean())*1e6,sub=(3,2,2),title=f"after shear compensation, rms {Z['Tot'].std():.2e}, max {J2['max_abs_tot']:.1e}",unit='ΔT/T ×1e6',fig=fig.number)
a=fig.add_subplot(3,2,3);L=np.arange(1,41);a.semilogy(L,J2['amp_bg'][1:41],'o-',label='background');a.semilogy(L,np.maximum(J2['amp_tot'][1:41],1e-12),'s-',label='background + bodies');a.axhline(1e-5,color='r',ls='--',label='1e-5');a.set_xlabel('ℓ');a.set_ylabel('rms ΔT/T per ℓ');a.legend(fontsize=7)
a=fig.add_subplot(3,2,4);k=np.abs(M)>0.01*np.abs(M).max();a.scatter(y[k&(M>0),0],y[k&(M>0),1],s=np.abs(M[k&(M>0)])/np.abs(M).max()*20,c='r',label=f"clusters {J2['n_clusters']}");a.scatter(y[k&(M<0),0],y[k&(M<0),1],s=np.abs(M[k&(M<0)])/np.abs(M).max()*20,c='b',label=f"voids {J2['n_voids']}");a.set_aspect('equal');a.set_xlabel('x (Gpc)');a.set_ylabel('y (Gpc)');a.set_title('bodies (projected), size ∝ |M|');a.legend(fontsize=7)
a=fig.add_subplot(3,2,5);Pf,Lt,D=Z3['Pf'],Z3['Lt'],Z3['D'];ll=np.arange(len(D));a.errorbar(Pf[:29,0],Pf[:29,1],[Pf[:29,2],Pf[:29,3]],fmt='o',color='grey',ms=3,label='Planck');a.plot(ll[2:31],Lt[2:31],'k',label='ΛCDM')
f2=J4['l2only']['best_f'];Dc=D.copy();Dc[2]=Z3['Df'][2]*f2**2+Z3['Dfix_new'][2];a.plot(ll[2:31],D[2:31],'r',label='CoE field + compensated bg');a.plot(ll[2:31],Dc[2:31],'m--',label=f'+ bodies cancel field ℓ=2 (f={f2})');a.set_xlabel('ℓ');a.set_ylabel('D_ℓ (μK²)');a.legend(fontsize=7)
a=fig.add_subplot(3,2,6);a.errorbar(l,yA,s,fmt='ko',ms=2);a.plot(l,m0,'g',label=f"η-edge locked χ²={mig['chi2_before']:.3f} r²={mig['r2_before']:.4f}");a.plot(l,m1,'m--',label=f"+ bodies χ²={mig['chi2_after']:.3f} r²={mig['r2_after']:.4f}");a.set_xlabel('l (deg)');a.set_ylabel('A/A_all');a.legend(fontsize=7)
fig.subplots_adjust(hspace=0.45,wspace=0.3,left=0.07,right=0.97,top=0.96,bottom=0.06);fig.savefig('s5_joint_summary.png',dpi=110)
