import json,sys,numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
sys.path.insert(0,'/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/rotfit');import fit
T0=2.7255e6;B='/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/'
web=np.load(B+'reproduce_cmb/webrec_D42.npy');ns=hp.npix2nside(len(web));V=np.array(hp.pix2vec(ns,np.arange(len(web)))).T
o=np.load(B+'quad_balance/mass_family_scan/runs/soft_2_30.npz')['o'];o=o[np.argsort(o[:,0])]
psi=np.degrees(np.arccos(np.clip(V@fit.vec(0,15),-1,1)));sm=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));sm=sm/sm.mean()-1
sm=np.asarray(hp.remove_dipole(sm));web=np.asarray(hp.remove_dipole(web-web.mean()));tot=sm+web
LMAX=600;l=np.arange(LMAX+1)
Dl=lambda m:l*(l+1)*hp.anafast(m,lmax=LMAX)*T0**2/(2*np.pi)
Dsm,Dw,Dt=Dl(sm),Dl(web),Dl(tot)
P=np.loadtxt(B+'reproduce_cmb/planckTT.txt');pl=P[:,0].astype(int);PD=P[:,1]
# check vs earlier per-l rms numbers: A_l=sqrt((2l+1)C_l/4pi)
Cl=hp.anafast(sm,lmax=6);print('check smooth per-l rms (1e-5):',np.round(np.sqrt((2*np.arange(7)+1)*Cl/(4*np.pi))/1e-5,2))
print('web map rms uK',np.std(web)*T0,' total rms uK',np.std(tot)*T0)
pd=dict(zip(pl,PD));rows=[]
print(f"{'l':>4} {'model D_l':>12} {'smooth':>12} {'web':>10} {'Planck':>9} {'ratio':>9}")
for L in list(range(2,13))+[20,30,50,100,200,300,500]:
    rows.append(dict(l=L,model=float(Dt[L]),smooth=float(Dsm[L]),web=float(Dw[L]),planck=float(pd.get(L,np.nan))))
    print(f"{L:4d} {Dt[L]:12.4g} {Dsm[L]:12.4g} {Dw[L]:10.4g} {pd.get(L,np.nan):9.1f} {Dt[L]/pd.get(L,np.nan):9.3g}")
json.dump(dict(rows=rows,lmax=LMAX,nside=ns),open('power_spectrum.json','w'),indent=1)
np.savetxt('model_Dl.txt',np.c_[l,Dt,Dsm,Dw],header='l D_l_total D_l_smooth D_l_web  [uK^2], D_l=l(l+1)C_l/2pi, C_l=anafast of dT/T map * T0^2')
fig,ax=plt.subplots(1,2,figsize=(18,7))
for a,xs in zip(ax,('log','linear')):
    a.plot(l[2:],Dt[2:],'k',lw=1.2,label='our model: smooth (soft p=2, M=30) + web, dipole removed');a.plot(l[2:],Dsm[2:],'C3--',lw=1,label='smooth ray-traced sky only');a.plot(l[2:],Dw[2:],'C0:',lw=1,label='recombination-era web only')
    a.errorbar(pl,PD,yerr=[P[:,2],P[:,3]],fmt='.',c='gold',ms=3,elinewidth=.4,label='Planck 2018 TT')
    a.set_yscale('log');a.set_xscale(xs);a.set_xlim(2,LMAX);a.set_ylim(1e-1,1e9);a.set_xlabel('multipole l');a.set_ylabel('D_l = l(l+1)C_l/2π  [μK²]');a.legend(fontsize=9)
ax[0].annotate('quadrupole',(2,Dt[2]),(4,Dt[2]*3),arrowprops=dict(arrowstyle='->'))
fig.suptitle('Toy-universe CMB angular power spectrum vs Planck 2018');plt.savefig('power_spectrum.png',dpi=110,bbox_inches='tight')
