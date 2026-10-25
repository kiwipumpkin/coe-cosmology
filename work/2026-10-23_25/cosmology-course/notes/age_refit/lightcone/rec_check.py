"""Geometry check, Chris's fixed-R_rec construction (illustrative, flat space, CoE frame, c=1).
ASSUMPTIONS: CoE at origin, static; R_rec=1 fixed. Every element moves radially; after crossing R_rec it follows the SAME law
r(s), s=time since its crossing: dr/ds=U, dU/ds=-mu/r^2 (Newtonian central mass, mu=G M/c^2 in R_rec units), U(0)=U_rec>1.
Elements cross R_rec continuously (earlier crossers are farther out). Each element emits once, at 3000 K, at its crossing event.
Inside R_rec = plasma, opaque. Photons travel straight at c in the CoE frame. We crossed at t=0 on the +x axis; now t=t*, D=r(t*).
For a sky direction theta (0 = away from CoE): the past ray from us must meet a crossing event; the first point it meets on the R_rec
circle (from outside) is that event, emitted at time t_e = t* - (path length). Rays that miss the R_rec disk meet no crossing event."""
import numpy as np, json
from scipy.integrate import solve_ivp
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
Rrec,Urec,mu,ts=1.0,3.0,1.0,2.0
sol=solve_ivp(lambda s,y:[y[1],-mu/y[0]**2],(0,ts),[Rrec,Urec],dense_output=True,rtol=1e-11,atol=1e-12)
r=lambda s: sol.sol(s)[0]; U=lambda s: sol.sol(s)[1]; D=r(ts); Rv=ts
rows=[]
for th in np.arange(0,360,30.):
    n=np.array([np.cos(np.radians(th)),np.sin(np.radians(th))])     # direction we look; photon came from +n side
    P=np.array([D,0.]); b=P@n; c=P@P-Rrec**2; disc=b*b-c
    row=dict(theta=th)
    if disc>=0 and -b-np.sqrt(disc)>0:
        L=-b-np.sqrt(disc); E=P+L*n; te=ts-L; alpha=np.degrees(np.arctan2(E[1],E[0]))%360
        row.update(hit=True,path=L,t_e=te,emit_xy=E.tolist(),alpha_CoE=alpha,beyond_CoE=bool(E[0]<0),elem_now_r=float(r(ts-te)) if ts-te<=ts else None,
                   photon_crosses_plasma=False)
    else: row.update(hit=False)
    rows.append(row)
cap=np.degrees(np.arcsin(Rrec/D))
fig,ax=plt.subplots(figsize=(19,12)); ax.set_aspect('equal'); c_=np.linspace(0,2*np.pi,400)
ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe0e0')); ax.plot(Rrec*np.cos(c_),Rrec*np.sin(c_),'r--',lw=3,label='R_rec (fixed); inside = plasma')
ax.plot(D+Rv*np.cos(c_),Rv*np.sin(c_),c='orange',lw=3,label=f'our visible circle R=c t*={Rv}')
for rr_,l in ((D-Rv,'D−R'),(D+Rv,'D+R')): ax.plot(rr_*np.cos(c_),rr_*np.sin(c_),':',c='gray',lw=2); ax.text(rr_*0.71,rr_*0.71,l,fontsize=16)
ax.plot([Rrec,D],[0,0],'b-',lw=4,label=f'our path since crossing (U {Urec}c → {U(ts):.2f}c)'); ax.plot(D,0,'b*',ms=26,label=f'us now, D={D:.2f}')
for q in rows:
    if q['hit']:
        E=q['emit_xy']; ax.annotate('',xy=(D,0),xytext=E,arrowprops=dict(arrowstyle='->',color='C2',lw=3)); ax.plot(*E,'o',c='C2',ms=12)
        ax.annotate(f"θ={q['theta']:.0f}°: path {q['path']:.2f} > c t*={ts}\nt_e={q['t_e']:.2f} < 0 (before our crossing): outside visible circle",E,fontsize=14,xytext=(-200,10 if E[1]>=0 else -25),textcoords='offset points')
    else:
        n=np.array([np.cos(np.radians(q['theta'])),np.sin(np.radians(q['theta']))]); ax.plot([D,D+Rv*n[0]],[0,Rv*n[1]],'k:',lw=1.5)
        ax.text(D+Rv*1.05*n[0],Rv*1.05*n[1],f"{q['theta']:.0f}°: no event",fontsize=13,ha='center')
ax.plot(0,0,'k*',ms=22); ax.text(0.05,-0.25,'CoE',fontsize=18)
nh=sum(q['hit'] for q in rows); nb=sum(q.get('beyond_CoE',False) for q in rows)
ax.set_title(f'Fixed R_rec, emission at crossing, straight photons (flat, CoE frame). Directions with an emission event: {nh}/12, beyond CoE: {nb}\n'
             f'Gap D−R_rec={D-Rrec:.2f} > c t*={ts}: no crossing event inside our visible circle (we outrun light, U>c). Disc of R_rec seen from us: asin(R_rec/D) of angular radius asin(R_rec/D) = {cap:.1f}°',fontsize=17)
ax.set_xlim(-D-Rv-0.3,D+Rv+0.8); ax.tick_params(labelsize=16); ax.grid(); ax.legend(fontsize=14,loc='lower left'); ax.set_xlabel('x / R_rec',fontsize=18)
fig.tight_layout(); fig.savefig('rec_check.png',dpi=80)
json.dump(dict(Rrec=Rrec,Urec=Urec,mu=mu,tstar=ts,D=D,U_now=U(ts),cap_deg=cap,rows=rows),open('rec_check.json','w'),indent=1,default=float)
print(D,U(ts),cap); [print(q) for q in rows]
