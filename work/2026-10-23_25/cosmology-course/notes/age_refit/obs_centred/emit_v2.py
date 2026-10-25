import pickle,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
d=pickle.load(open('oc_m253_rays.pkl','rb')); J=json.load(open('oc_m253.json')); S=J['sync_chris']; RD=d['RD']
rows={r['ell']:r for r in S['rows']}; vis=[];cut=[];stop=[]
for e,rd in zip(d['ells'],RD):
    r=rows.get(float(e))
    if r is None or 'T_obs_K' not in r: stop.append(e); continue
    p=np.radians(r['phi_e_deg']); x=r['R_e_over_D']*np.cos(p)-1; y=r['R_e_over_D']*np.sin(p)
    if x< -1 or int(rd['nc'].max())>0: cut.append((e,x,y))
    else: vis.append((e,x,y,r['T_obs_K']))
fig,a=plt.subplots(figsize=(9,7)); v=np.array(vis)
sc=a.scatter(v[:,1],v[:,2],c=v[:,3],cmap='coolwarm',s=40,zorder=3); plt.colorbar(sc,label='observed T (K)')
for e in [c[0] for c in cut]+stop:
    th=np.radians(e); col='k' if e in [c[0] for c in cut] else 'm'
    a.plot([0,0.25*np.cos(th)],[0,0.25*np.sin(th)],c=col,lw=2 if col=='k' else 1,ls='-' if col=='k' else '--')
a.plot([],[],'k-',label=f'horizon-cut direction ({len(cut)})'); a.plot([],[],'m--',label=f'integrator-stopped, no emitter ({len(stop)})')
a.plot(0,0,'b*',ms=14,label='us'); a.plot(-1,0,'kx',ms=12,label='CoE'); a.axvline(-1,c='grey',lw=.6,ls=':')
a.set_aspect('equal'); a.set_xlabel('x / D'); a.set_ylabel('y / D'); a.legend(loc='upper right',fontsize=8)
a.set_title(f"Where today's CMB was emitted (M=2.5e53 kg, η fwd, U_o={J['Uo_c']:.2f}c)\nvisible {len(vis)} / cut {len(cut)} / stopped {len(stop)} of {len(RD)} directions (6°)")
fig.text(.02,.01,"Rules: emission only at t≥t* (our 3000 K crossing, t*≈24.8 Myr) from matter recombined by then; exact-GR rays; T_obs=T_emit/(1+z), mean 2.725 K.\n"
 "CoE side is a horizon: emitters at x<-1 or rays that pass through the centre are cut (no emitter). Integrator-stopped rays get no point.",fontsize=7.5)
plt.tight_layout(rect=(0,.06,1,1)); plt.savefig('emit_v2.png',dpi=90)
print(len(vis),len(cut),len(stop),[c[0] for c in cut],stop, J['Uo_c'])
