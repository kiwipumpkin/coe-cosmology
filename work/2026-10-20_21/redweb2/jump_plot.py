import json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('plot4.py').read());J=json.load(open('thin4_jump.json'));S=json.load(open('thin4_jump_step.json'));T=json.load(open('thin4_results.json'))
c=[o for o in S if o['gid']==21717][0];z=np.load('jump_data.npz');k=np.where(z['gid']==21717)[0];print('u',z['u'][k],'M',c['M_parts'],'l,b',np.degrees(np.arctan2(c['pos'][1],c['pos'][0]))%360,np.degrees(np.arcsin(c['pos'][2]/c['r'])))
R=J['rows'];f=[r['f']*100 for r in R]
fig,a=plt.subplots(2,1,figsize=(10,8),sharex=True);a[0].plot(f,[r['chi2'] for r in R],'ko-');a[0].set_ylabel('χ² vs Planck (ℓ=2–30)');a[0].axvline(1.705,color='r',ls=':',label='jump: void 21717 added (u=0.01700…)');a[0].legend()
for L,cc in((20,'tab:blue'),(25,'tab:green'),(30,'tab:orange')):a[1].plot(f,[r['Dl'][L] for r in R],'o-',color=cc,label=f'D_ℓ at ℓ={L}')
a[1].set_ylabel('D_ℓ (μK²)');a[1].set_xlabel('kept fraction (%) — draw 0, 4 Gpc read');a[1].legend();a[1].axvline(1.705,color='r',ls=':');plt.tight_layout();plt.savefig('thin4_jump_chi2.png',dpi=110);plt.close()
plt.figure(figsize=(11,7));l=np.arange(2,31);B={r['f']:r for r in R}
for lab,D,cc,m in [('1.5% (draw 0)',B[0.015]['Dl'],'tab:blue','o'),('1.70% just before',B[0.017]['Dl'],'tab:green','s'),('1.71% just after',B[0.0171]['Dl'],'tab:red','^'),('1.75% (draw 0)',B[0.0175]['Dl'],'tab:purple','D')]:
  D=np.array(D);plt.plot(l,D[2:31],color=cc,marker=m,ms=5,label=f'{lab}: χ²={chiT(D) if False else [r for r in R if r["Dl"]==list(D)][0]["chi2"]:.1f}, r²={r2(D):.2f}')
plt.errorbar(Pl[:,0],Pl[:,1],yerr=[Pl[:,2],Pl[:,3]],fmt='o',color='k',ms=4,capsize=2,label='Planck TT');plt.ylim(0,3000);plt.xlabel('ℓ');plt.ylabel('D_ℓ (μK²)');plt.title('Jump between 1.70% and 1.71% kept (4 Gpc read)');plt.legend();plt.tight_layout();plt.savefig('thin4_jump_spectra.png',dpi=110);plt.close()
fig,ax=plt.subplots(1,2,figsize=(16,7.5));p=np.array(c['pos'])
for a_,(i,j,lab) in zip(ax,[(0,1,'x–y'),(0,2,'x–z')]):
  t=np.linspace(0,2*np.pi,400);a_.plot(4*np.cos(t),4*np.sin(t),'--',color='tab:purple',label='4 Gpc sampling sphere (great circle)')
  rs=np.sqrt(max(16-p[3-i-j]**2,0));a_.plot(rs*np.cos(t),rs*np.sin(t),':',color='tab:purple',label=f'4 Gpc sphere slice at body’s {"xyz"[3-i-j]}={p[3-i-j]:.2f}')
  X=z['X'];s=(z['u']>=0.0170)&(z['u']<0.0171);o=~s
  a_.scatter(X[o,i],X[o,j],s=4,c='0.7',label='other bodies added 1.5→1.75%');a_.scatter(X[s,i],X[s,j],s=15,c='tab:blue',label='bodies added at 1.70→1.71% step')
  a_.scatter(p[i],p[j],s=200,marker='*',c='red',zorder=5,label=f'void 21717 (r={c["r"]:.3f} Gpc, {1000*c["dsphere"]:.0f} Mpc outside sphere)');a_.plot(0,0,'k+',ms=14,mew=2)
  a_.set_xlim(-10,10);a_.set_ylim(-10,10);a_.set_aspect('equal');a_.set_xlabel(f'galactic {"xyz"[i]} (Gpc)');a_.set_ylabel(f'galactic {"xyz"[j]} (Gpc)');a_.set_title(lab+' projection')
ax[0].legend(fontsize=8,loc='lower right');plt.tight_layout();plt.savefig('thin4_jump_bodies.png',dpi=110)
