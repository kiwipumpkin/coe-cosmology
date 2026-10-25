import sys;k=int(sys.argv[1]);sys.argv=['x'];src=open('flowR.py').read().split('\nif __name__')[0];exec(src)
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
ns=64;V=np.array(hp.pix2vec(ns,np.arange(hp.nside2npix(ns)))).T
Rf,Rq=grids(500000+k);Rq=Rq+tocoarse(Rf);fq=phi(Gq,Rq,Gq['dx']);ff_=phi(Gf,Rf,Gf['dx']);fs=phi(Gf,Rf,Gq['dx']);ff=[a-b for a,b in zip(ff_,fs)];del ff_,fs
T0=2.7255e6;LM=110;Ll=np.arange(LM+1);fac=Ll*(Ll+1)*T0**2/(2*np.pi)
sm=hp.alm2map(hp.map2alm(np.load('/workspace/redweb2/sm_red.npy'),lmax=95),ns,lmax=95)
th,ph=hp.pix2ang(ns,np.arange(hp.nside2npix(ns)));et=np.stack([np.cos(th)*np.cos(ph),np.cos(th)*np.sin(ph),-np.sin(th)],1);ep=np.stack([-np.sin(ph),np.cos(ph),0*ph],1)
g=lambda m:hp.alm2map_der1(hp.map2alm(m,lmax=191),ns,lmax=191)[1:];gs=g(sm)
P=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-full_R3.01.txt');Pb=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-TT-binned_R3.01.txt');Pb=Pb[Pb[:,0]<=LM+5]
Th=np.loadtxt('/workspace/planck/COM_PowerSpect_CMB-base-plikHM-TTTEEE-lowl-lowE-lensing-minimum-theory_R3.01.txt');LT=dict(zip(Th[:,0].astype(int),Th[:,1]))
m30=(P[:,0]>=2)&(P[:,0]<=30);lp=P[m30,0].astype(int);dp_=P[m30,1]
def st(D):x=D[lp];s=np.where(x>dp_,P[m30,3],P[m30,2]);return float(np.sum(((x-dp_)/s)**2)),float(np.corrcoef(x,dp_)[0,1]**2)
def Dl(m):return fac*hp.anafast(np.asarray(hp.remove_dipole(m)),lmax=LM)
E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];pm=np.median(P[(P[:,0]>=5)&(P[:,0]<=30),1])
Lc=np.array([0,0]+[LT[q] for q in range(2,LM+1)]);D0=Dl(sm)
fig,ax=plt.subplots(4,1,figsize=(12,22))
for a,R in zip(ax,(2,4,8,14)):
    o=comp(fq,ff,float(R));gw=g(o['sw_grow']);al0=o['al_grow']
    cq,r=st(Lc);a.plot(Ll[2:],Lc[2:],'k-',lw=2,label=f'ΛCDM χ²={cq:.1f}, r²={r:.2f}')
    cq,r=st(D0);a.plot(Ll[2:],D0[2:],'r--',lw=2,label=f'smooth red χ²={cq:.1f}, r²={r:.2f}')
    med=[]
    for j,e in enumerate(E):
        dt=gs[0]+e*gw[0];dp=gs[1]+e*gw[1];al=e*al0+o['ac']
        mp=sm+e*o['sw_grow']+o['rs_cc']+e*(o['isw_grow']+o['mlf_grow'])+(al*et).sum(1)*dt+(al*ep).sum(1)*dp
        D=Dl(mp);cq,r=st(D);med.append(np.median(D[5:31]));a.plot(Ll[2:],D[2:],color=plt.cm.viridis(j/(len(E)-1)),label=f'E={e:g} χ²={cq:.4g}, r²={r:.2f}')
    med=np.array(med);i=np.where((med[:-1]>=pm)&(med[1:]<pm))[0]
    br=f'E {E[i[0]+1]:g}–{E[i[0]]:g} (medians {med[i[0]+1]:.0f}, {med[i[0]]:.0f})' if len(i) else 'none'
    print(f'u{k} R={R}: bracket {br}; medians',np.round(med,1).tolist(),flush=True)
    lo=P[P[:,0]<30];a.errorbar(lo[:,0],lo[:,1],yerr=[lo[:,2],lo[:,3]],fmt='o',color='k',ms=4,label='Planck')
    a.errorbar(Pb[:,0],Pb[:,1],yerr=[Pb[:,2],Pb[:,3]],fmt='s',color='k',mfc='w',ms=4,label='Planck binned (ℓ≥30)')
    a.axvspan(100,LM,color='0.85',alpha=.6);a.set_xlim(2,LM);a.set_yscale('log');a.set_ylim(1e-3,None);a.set_ylabel('D_ℓ (μK²)')
    a.set_title(f'web riding CoE flow, grow_full, R={R} Gpc, #{k} (nside 64): Planck-level bracket {br}');a.legend(loc='lower center',ncol=3,fontsize=7,framealpha=.85)
ax[-1].set_xlabel('ℓ');fig.tight_layout();fig.savefig(O+f'ladder_full_u{k}.png',dpi=90)
