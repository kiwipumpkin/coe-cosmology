import json,numpy as np,sys;sys.path.insert(0,'.')
exec(open('rec_fit6.py').read().split('BOX=')[0])
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':20})
o=json.load(open('rec_fit6.json'));K=list(o)
# geometry
P=o['A_TD']['params'];x=o['A_TD']['x'];U0=P['U0'];mu=P['mu'];ts=P['ts'];sol=traj(U0,mu,ts)
fig,ax=plt.subplots(1,2,figsize=(26,11));tt=np.linspace(0,ts,300)
ax[0].plot(tt,sol(tt)[0],lw=4,label='D(t) (our distance from CoE)');ax[0].plot(tt,tt,lw=4,label='R = c t* (visible radius)')
ax[0].plot(tt,sol(tt)[1]*100,'--',lw=3,label='U(t)×100 [c]');ax[0].set_xlabel('t since we crossed R_rec');ax[0].legend();ax[0].grid()
c_=np.linspace(0,2*np.pi,300);a=ax[1];a.set_aspect('equal')
for tf,cl in ((ts/3,'C0'),(2*ts/3,'C1'),(ts,'C3')):
    d,R=sol(tf)[0],tf;a.plot(d+R*np.cos(c_),R*np.sin(c_),c=cl,lw=3,label=f't*={tf:.0f}: D={d:.0f}, R={R:.0f}')
    for rr in (d-R,d+R):a.plot(rr*np.cos(c_),rr*np.sin(c_),':',c=cl,lw=2)
    a.plot(d,0,'*',c=cl,ms=22)
a.plot(0,0,'k*',ms=24);a.text(10,-90,'CoE',fontsize=22);a.legend(fontsize=16,loc='upper left');a.set_xlim(-1900,1900);a.set_ylim(-1000,1000)
a.set_title('Visible circle R=ct* (solid) & CoE shells D±R (dotted)')
fig.suptitle(f'rec_fit6 geometry (A, with TD best fit): R grows at c, D grows at U>c (U0={U0:.1f}c→{P["U"]:.2f}c)',fontsize=24)
fig.tight_layout();fig.savefig('rec_fit6_geometry.png',dpi=60)
# T(psi)
fig,ax=plt.subplots(2,2,figsize=(26,16))
for a,k in zip(ax.flat,K):
    r=o[k]['1deg'];e=np.arange(360);a.plot(e,np.array(r['f'])*100,lw=3,label='T/⟨T⟩−1 raw');a.plot(e,np.array(r['r'])*100,lw=3,label='dipole removed')
    a.set_title(f"{k}: ⟨T⟩={r['Tmean']:.3f} K, dip={r['dipole']*100:.2f}%, quad={r['quad']*100:.2f}%",fontsize=20);a.set_xlabel('ψ from CoE direction [deg]');a.set_ylabel('%');a.grid();a.legend()
fig.tight_layout();fig.savefig('rec_fit6_T.png',dpi=60)
# Migkas
ang=np.arange(360.);fig,ax=plt.subplots(2,2,figsize=(26,16))
for a,k in zip(ax.flat,K):
    r=o[k]['1deg'];sh=F.shifted(np.array(r['af']),r['l0']);a.plot(ang,sh,lw=3,label='model')
    a.errorbar(np.arange(360)[F.IMIG] if np.ndim(F.IMIG) else ang,F.Y,yerr=1/np.sqrt(F.W),fmt='ko',ms=5,label='Migkas A/A_all')
    a.set_title(f"{k}: χ²={r['chi2']:.1f}, σ={r['sigma']:.2f}, ℓ0={r['l0']:.0f}°",fontsize=20);a.set_xlabel('ℓ [deg]');a.grid();a.legend()
fig.tight_layout();fig.savefig('rec_fit6_migkas.png',dpi=60)
# summary table
fig,a=plt.subplots(figsize=(30,12));a.axis('off');rows=[]
for k in K:
    p=o[k]['params'];r6,r1=o[k]['6deg'],o[k]['1deg']
    rows.append([k,f"{p['U0']:.2f}",f"{p['mu']:.3g}",f"{p['ts']:.0f}",f"{p['D']:.0f}",f"{p['U']:.2f}",f"{p['2mu_over_min_re']:.2f}",f"{r6['chi2']:.1f}/{r6['sigma']:.2f}",f"{r1['chi2']:.1f}/{r1['sigma']:.2f}",f"{r1['Tmean']:.3f}",f"{r1['dipole']*100:.2f}",f"{r1['quad']*100:.2f}",f"{r1['rms_nodip']*100:.3f}",
      ' '.join(f"{t}[{v[0]:.3f},{v[2]:.3f}]" for t,v in r1['terms'].items() if t in('grav','gamma','dop'))+f" str={r1['terms']['stretch'][1]:.0f}"])
t=a.table(cellText=rows,colLabels=['run','U0','μ','t*','D','U now','2μ/r_min','Mig χ²/σ 6°','Mig χ²/σ 1°','⟨T⟩ K','dip %','quad %','rms nodip %','terms (min,max)'],loc='center');t.auto_set_font_size(False);t.set_fontsize(15);t.scale(1,4);t.auto_set_column_width(range(14))
a.set_title('rec_fit6: A=wavelength-only stretch, B=no stretch; TD = grav clock + γ',fontsize=24);fig.savefig('rec_fit6_summary.png',dpi=60,bbox_inches='tight')
