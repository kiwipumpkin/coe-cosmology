import numpy as np,healpy as hp,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
O='/workspace/redweb2/realweb/los_flow_R/';T0=2.7255e6;E=0.0131
P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt');Pb=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-binned_R3.01.txt')
Th=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt');LT=dict(zip(Th[:,0].astype(int),Th[:,1]))
def sky(z,ns):
    lm=3*ns-1;th,ph=hp.pix2ang(ns,np.arange(hp.nside2npix(ns)))
    et=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1);ep=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
    g=lambda m:hp.alm2map_der1(hp.map2alm(m,lmax=lm),ns,lmax=lm)[1:]
    gs=g(z['sm']);gw=g(z['sw_grow']);dt=gs[0]+E*gw[0];dp=gs[1]+E*gw[1];al=E*z['al_grow']+z['ac']
    out=z['sm']+E*z['sw_grow']+z['rs_cc']+E*(z['isw_grow']+z['mlf_grow'])
    return out+(al*et).sum(1)*dt+(al*ep).sum(1)*dp
S={}
for ns in (32,64,128,256):
    try:z=dict(np.load(O+f'hl_R4_ns{ns}.npz'))
    except Exception:continue
    lm=3*ns-1;c=hp.anafast(np.asarray(hp.remove_dipole(sky(z,ns))),lmax=lm);L=np.arange(lm+1);S[ns]=L*(L+1)*c*T0**2/(2*np.pi)
for ns in S:print(ns,'D_l at l=5,10,20,40,60,80,100,150,200:',[round(float(S[ns][q]),1) if q<len(S[ns]) else None for q in (5,10,20,40,60,80,100,150,200)])
ks=sorted(S)
for a,b in zip(ks[:-1],ks[1:]):
    r=np.abs(S[a][2:len(S[a])]/S[b][2:len(S[a])]-1);bad=np.where(r>0.1)[0];print(f'ns{a} vs ns{b}: first l with >10% diff:',bad[0]+2 if len(bad) else None)
np.save(O+'hl_spectra.npy',{k:v for k,v in S.items()},allow_pickle=True)
D=S[128];lm=len(D)-1;Ls=np.arange(2,lm+1)
print('D_l 30..120 step10',[round(float(D[q]),1) for q in range(30,121,10)],'peak',D[2:].max())
def stats(model,lo,hi):
    m=(P[:,0]>=lo)&(P[:,0]<=hi);l_=P[m,0].astype(int);d=P[m,1];x=model(l_);s=np.where(x>d,P[m,3],P[m,2])
    return float(np.sum(((x-d)/s)**2)),len(l_),float(np.corrcoef(x,d)[0,1]**2)
ours=lambda l_:D[l_];lc=lambda l_:np.array([LT[q] for q in l_])
for lo,hi in ((2,30),(2,107),(2,383)):
    for nm,f in (('CoE #113 R=4 E=0.0131',ours),('LCDM Planck18 best fit',lc)):
        c2,n,r=stats(f,lo,hi);print(f'l {lo}-{hi}: {nm}: chi2={c2:.1f} (N={n}), Pearson r2={r:.3f}')
c2a,_,ra=stats(ours,2,30);c2l,_,rl=stats(lc,2,30)
fig,ax=plt.subplots(2,1,figsize=(12,10))
for a,(x0,x1,log) in zip(ax,((2,2508,True),(2,200,False))):
    L2=np.arange(2,2509);a.plot(L2,[LT[q] for q in L2],'k-',lw=2,label=f'ΛCDM (Planck 2018 best fit) χ²(ℓ2–30)={c2l:.1f}, r²={rl:.2f}')
    a.plot(Ls,D[2:],'-',color=plt.cm.viridis(.3),lw=2,label=f'CoE #113 grow_full R=4 E=0.0131 (nside128) χ²(ℓ2–30)={c2a:.1f}, r²={ra:.2f}')
    a.axvspan(107,383,color='0.85',alpha=.5,label='ours beyond nside-converged range (ℓ>107)')
    lo=P[P[:,0]<30];a.errorbar(lo[:,0],lo[:,1],yerr=[lo[:,2],lo[:,3]],fmt='o',color='k',ms=3,label='Planck 2018 TT (ℓ<30)')
    a.errorbar(Pb[:,0],Pb[:,1],yerr=[Pb[:,2],Pb[:,3]],fmt='s',color='tab:red',ms=3,label='Planck 2018 TT binned (ℓ≥30)')
    if log:a.set_xscale('log')
    a.set_xlim(x0,x1);a.set_ylim(0,6500 if log else 3000);a.set_ylabel('D_ℓ (μK²)');a.set_xlabel('ℓ');a.legend(loc='upper left' if log else 'upper right',fontsize=8,framealpha=.85)
ax[0].set_title('TT power spectrum: CoE model vs ΛCDM vs Planck (full range, log ℓ)');ax[1].set_title('zoom ℓ 2–200 (linear)')
fig.tight_layout();fig.savefig(O+'highl_spectrum.png',dpi=90)
