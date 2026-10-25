import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
law,mode,tag=sys.argv[1],sys.argv[2],(sys.argv[3] if len(sys.argv)>3 else '')
if tag=='_gtc': os.environ['VMODE']='gtc'
os.environ.setdefault('GRZ','gr'); sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
x=np.array(json.load(open(f'fwd/lc_{law}_{mode}{tag}_ckpt.json'))['x']); P,o,af,chi,l0=evl(x)
L_=F.local(o); ell=np.arange(360.)
def nod(S):
    th=np.radians(ell); X=np.c_[np.ones(360),np.cos(th),np.sin(th)]; c=np.linalg.lstsq(X,np.log(S),rcond=None)[0]; return S/np.exp(c[1]*np.cos(th)+c[2]*np.sin(th))
def orient(S): return np.interp((ell-l0)%360,np.arange(361.),np.append(S,S[0]))
fig,ax=plt.subplots(2,1,figsize=(16,10))
sigma=0.5*(F.MLO+F.MHI); A=F.shifted(af,l0)
ax[0].errorbar(F.IMIG,F.Y,yerr=sigma,fmt='o',ms=6,c='k',label='Migkas 71 bins')
ax[0].plot(ell,A,lw=3,c='C3',label=f'model χ²={chi:.2f} (null 45.1), l0={l0:.0f}°')
ax[0].set_xlim(0,360);ax[0].set_ylabel('A/A_all',fontsize=18);ax[0].legend(fontsize=14)
cs=[]
c=F.cmb_sr_local(o,0.23); S1=c['S']; cu=F.cmb_uniform(o,1/1100.); S2=cu['S']
for S,nm,col in ((S1,'local-SR q=0.23','C0'),(S2,'uniform surface a=1/1100','C2')):
    if not np.all(np.isfinite(S)): cs.append(f'{nm}: incomplete sky ({int(np.isfinite(S).sum())}/360)'); continue
    So=orient(S); ax[1].plot(ell,So/So.mean()-1,lw=3,c=col,label=nm+' dT/T')
    Sn=orient(nod(S)); ax[1].plot(ell,Sn/Sn.mean()-1,lw=2,ls='--',c=col,label=nm+' m=1 removed')
    r=CC.resid(CC.G_FULL,So); ax[1].plot(ell,r,lw=2,ls=':',c=col,label=nm+f' G·S resid ({100*(1-CC.rms(r)/S0):.1f}% cancel)')
ax[1].plot(ell,CC.G_FULL-1,c='k',lw=1.5,label='G−1 (Migkas curve)')
ax[1].set_xlim(0,360);ax[1].set_xlabel('galactic longitude ℓ (deg)',fontsize=18);ax[1].set_ylabel('dT/T, residual',fontsize=18);ax[1].legend(fontsize=11,ncol=2)
for a in ax: a.tick_params(labelsize=14);a.grid()
ttl=f"{law} {mode}{tag}  U_o/c={L_['Uo']/F.C:.3g}  grav/kin={L_['grav_over_kin']:.3g}"+('  | '+'; '.join(cs) if cs else '')
fig.suptitle(ttl,fontsize=17);fig.tight_layout();out=f'fwd/curves_{law}_{mode}{tag}.png';fig.savefig(out,dpi=80);print(out)
