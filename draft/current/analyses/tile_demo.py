import json, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Circle
exec(open("scan2.py").read().split("out={")[0])
KG=7.7e51/0.0039; Rv=14.3
d=json.load(open("web12.json"))["12_1e+52"]
M=np.array(d["masses"]); A=np.radians(d["ang"]); r=np.array(d["dist_R"])*Rv
Dc,l0=58.2,np.radians(15.2)
obs=-Dc*np.array([np.cos(l0),np.sin(l0)])   # CoE at origin lies in direction l0 from observer
fx=obs[0]+r*np.cos(A); fy=obs[1]+r*np.sin(A)
# tile: alternating ±1e52 kg, neighbour spacing ~ mean fitted spacing
order=np.argsort(A); near=r<40
sp=np.median([np.hypot(fx[order][i]-fx[order][(i+1)%12],fy[order][i]-fy[order][(i+1)%12]) for i in range(12) if near[order][i] and near[order][(i+1)%12]])
rin,rout=Dc-Rv,Dc+Rv; tx=[];ty=[];tm=[]
rings=np.arange(rin+sp/2,rout,sp*0.87)
for j,rr in enumerate(rings):
    n=int(round(2*np.pi*rr/sp)); n+=n%2
    for i in range(n):
        a=2*np.pi*(i+0.5*(j%2))/n; x,y=rr*np.cos(a),rr*np.sin(a)
        if np.hypot(x-obs[0],y-obs[1])<34: continue
        tx.append(x);ty.append(y);tm.append(1e52*(-1)**(i+j))
tx,ty,tm=map(np.array,(tx,ty,tm))
# effect of the tiled extras on the fitted residual
base=R-sum(M[i]/KG*field(A[i],r[i]/Rv) for i in range(12))
dxo,dyo=tx-obs[0],ty-obs[1]; ta=np.arctan2(dyo,dxo); td=np.hypot(dxo,dyo)/Rv
extra=sum(tm[i]/KG*field(ta[i],td[i]) for i in range(len(tm)))
print("spacing Gpc",sp,"tiles",len(tm),"net tiles",tm.sum(),"rms fit only",rms(base),"with tiles",rms(base-extra))
fig,ax=plt.subplots(figsize=(9,9),dpi=140)
ax.add_patch(Circle((0,0),rout,color='#f1ece2',zorder=0)); ax.add_patch(Circle((0,0),rin,color='white',zorder=0))
s=lambda m:4.5*(abs(m)/1e52)**(1/3)+0.6
for x,y,m in zip(tx,ty,tm): ax.add_patch(Circle((x,y),s(m),color='#c0392b' if m>0 else '#2e86c1',alpha=0.35,lw=0,zorder=1))
for x,y,m in zip(fx,fy,M): ax.add_patch(Circle((x,y),s(m),color='#c0392b' if m>0 else '#2e86c1',alpha=0.95,ec='k',lw=1.2,zorder=3))
ax.add_patch(Circle(obs,Rv,fill=False,ls='--',lw=1.6,color='k',zorder=4))
ax.plot(0,0,'k*',ms=16,zorder=5); ax.text(2,-5,'CoE',fontsize=11)
ax.plot(*obs,'o',color='gold',mec='k',ms=10,zorder=5); ax.text(obs[0]-4,obs[1]-5,'observer',fontsize=10,ha='right')
ax.plot([0,obs[0]],[0,obs[1]],'k:',lw=1)
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([],[],marker='o',ls='',color='#c0392b',mec='k',ms=10,label='fitted cluster (≤1×10⁵² kg)'),
 Line2D([],[],marker='o',ls='',color='#2e86c1',mec='k',ms=10,label='fitted void (≥−1×10⁵² kg)'),
 Line2D([],[],marker='o',ls='',color='#c0392b',alpha=.35,ms=10,label='pattern repeated: cluster'),
 Line2D([],[],marker='o',ls='',color='#2e86c1',alpha=.35,ms=10,label='pattern repeated: void'),
 Line2D([],[],ls='--',color='k',label='vision sphere, 14.3 Gpc')],loc='upper right',fontsize=8.5,frameon=True)
L=max(rout,np.abs(np.r_[fx,fy]).max())+8; ax.set_xlim(-L,L); ax.set_ylim(-L,L); ax.set_aspect('equal'); ax.set_xlabel('Gpc'); ax.set_ylabel('Gpc')
ax.set_title(f'Fitted clusters and voids around our vision sphere (12 objects, residual {d["rms"]:.1e}, below 10⁻⁵),\nthen the alternating pattern repeated around the CoE shell ({rin:.1f}–{rout:.1f} Gpc), spacing ≈ {sp:.0f} Gpc',fontsize=10)
plt.tight_layout(); plt.savefig('tile_demo.png'); print('saved')
# refit the 12 near objects with the repeated pattern held fixed in the background
from scipy.optimize import least_squares
N=12; kc=1e52/KG; Rt=R-extra
def model(q):
    return sum(q[i]*field(np.radians(q[N+i]),q[2*N+i]) for i in range(N))
q0=np.r_[np.clip(M/KG,-kc*0.999,kc*0.999),np.degrees(A),np.clip(r/Rv,1.06,4.99)]; best=None
for s in range(3):
    rng=np.random.default_rng(s); qq=q0.copy()
    if s: qq[:N]*=rng.uniform(0.5,1.2,N); qq[N:2*N]+=rng.uniform(-8,8,N); qq[:N]=np.clip(qq[:N],-kc*0.999,kc*0.999)
    f=least_squares(lambda q:np.r_[model(q)-Rt,1e-3*np.sqrt(len(th))/(kc*N)*q[:N].sum()],qq,bounds=(np.r_[[-kc]*N,[-720]*N,[1.05]*N],np.r_[[kc]*N,[1080]*N,[5]*N]),max_nfev=600,x_scale=np.r_[[kc]*N,[30]*N,[0.5]*N])
    if best is None or f.cost<best.cost: best=f
q=best.x; print("refit with tiles fixed rms",rms(Rt-model(q)),"net",q[:N].sum()*KG,"at cap",int(np.sum(np.abs(q[:N])>0.98*kc)))
json.dump(dict(masses=(q[:N]*KG).tolist(),ang=(q[N:2*N]%360).tolist(),dist_R=q[2*N:].tolist(),rms=rms(Rt-model(q))),open("tile_refit.json","w"),indent=1)
for m,a,rr in sorted(zip(q[:N]*KG,q[N:2*N]%360,q[2*N:]*Rv),key=lambda t:t[1]): print(f"{a:6.1f} {rr:5.1f} Gpc {m:+.2e}")
