import numpy as np,pickle,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
P=pickle.load(open('xy.xy.pkl','rb'));o=np.load('xy.npz')['o']
import sys
Ro_e=0.0009162958051463858  # R(t_e, w=1.008) of psi=90 ray, proxy for observer shell
f,ax=plt.subplots(2,2,figsize=(13,11));a1,a2,a3,a4=ax.ravel();cm=plt.cm.viridis;rep=[]
for i,p in enumerate(P):
    c=cm(i/8);w,ph,t=p['w'],p['ph'],p['t'];R=p['R'];phs=ph[::20];ts=t[::20]
    a1.plot(w*np.cos(ph),w*np.sin(ph),color=c,label=f"ψ={p['psi']:.0f}°")
    a2.plot(R*np.cos(phs),R*np.sin(phs),color=c)
    sg=np.sign(-np.gradient(w));sg[sg==0]=1
    al=np.degrees(np.arccos(np.clip(sg*np.sqrt(p['c2']),-1,1)))  # angle propagation vs outward radial
    a3.plot(t/t[0],al,color=c)
    # comoving direction angle (geometric turning in w-phi plane, not metric)
    x,y=w*np.cos(ph),w*np.sin(ph);hd=np.degrees(np.unwrap(np.arctan2(-np.gradient(y),-np.gradient(x))))
    rep.append(dict(psi=float(p['psi']),alpha_obs=float(al[0]),alpha_emit=float(al[-1]),dalpha=float(al[0]-al[-1]),
      max_dalpha_step_t_frac=float(t[np.argmax(np.abs(np.gradient(al)))]/t[0]),
      half_turn_t_frac=float(t[np.argmin(np.abs(np.abs(al-al[-1])-0.5*abs(al[0]-al[-1])))]/t[0]),
      comoving_heading_turn_deg=float(hd[0]-hd[-1]),phi_emit_deg=float(np.degrees(ph[-1])),w_emit=float(w[-1])))
    if p['psi'] in (45,90,135):
        te=t[-1];Re=R[-1];Rn=np.interp(0,[0],[0]) 
        we,pe=w[-1],ph[-1]
        a4.plot([Re*np.cos(pe)],[Re*np.sin(pe)],'o',color=c,mfc='none');
        a4.plot([0,1.3e-3*np.cos(pe)],[0,1.3e-3*np.sin(pe)],'--',color=c,lw=.8);a4.plot(R*np.cos(phs),R*np.sin(phs),color=c,label=f"ψ={p['psi']:.0f}°")
        a4.plot([Re*np.cos(pe),Ro_e],[Re*np.sin(pe),0],':',color=c,lw=.8)
for a in (a1,a2,a4):
    a.plot(0,0,'k*',ms=14);a.text(0.002,0.0,'CoE',fontsize=10);a.plot(1,0,'ro',ms=7);a.grid(alpha=.3)
a1.set_aspect('equal');a1.set_ylim(-0.02,0.15);a1.set_xlim(-0.05,1.2)
a4.set_xlim(-0.0001,0.0013);a4.set_ylim(-0.00005,0.00015)
a4.plot(Ro_e*1,0,'rs',mfc='none',ms=7)
a1.set_title('(1) comoving paths: x=w cosφ, y=w sinφ (CoE left, observer red)');a1.legend(fontsize=8)
a2.set_title('(2) areal/physical paths: R(t,w) at each step\n(where the photon was at time t; R(t0,1)=1)')
a3.set_xlabel('t / t0');a3.set_ylabel('α [deg] (photon propagation vs outward radial)');a3.set_xscale('log');a3.grid(alpha=.3)
a3.set_title('(3) α along each ray (180° = moving straight toward CoE)')
a4.set_title('(4) zoom of (2) near emission: emitter at emission (open circle),\nred square = R at x-axis ≈ observer shell at t_e (uses ψ=90 emitter R, w=1.008), dashed = emitter radial line,\ndotted = straight chord emitter→observer-at-emission')
plt.tight_layout();plt.savefig('ray_paths.png',dpi=90)
json.dump(rep,open('ray_paths.json','w'),indent=1)
for r in rep:print({k:round(v,4) for k,v in r.items()})
print('R obs at emission',Ro_e, 'emit R', [p['R'][-1] for p in P][:3])
