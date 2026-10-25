import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
sys.path.insert(0,'/workspace/figs_letter_exact_python/makers');sys.path.insert(0,'/tmp')
from physics import rotate_xyz
from spots import SPOTS
M=json.load(open('/tmp/M25.json'));L=np.array(M['l'])
CR={'A':[22,99,262,310],'B':[14,104,266,321]}
def sgnfun(k):
    c=CR[k];mids=[(c[i]+((c[(i+1)%4]-c[i])%360)/2)%360 for i in range(4)]
    s=[np.sign(np.interp(m,L,np.array(M[k]),period=360)) for m in mids];return c,s
ZS={k:sgnfun(k) for k in 'AB'}
def zone(k,t):
    c,s=ZS[k];t=np.asarray(t)%360;out=np.full(t.shape,s[3])
    for i in range(4):
        lo,hi=c[i],c[(i+1)%4]
        m=((t-lo)%360)<((hi-lo)%360);out=np.where(m,s[i],out)
    return out
print({k:ZS[k] for k in 'AB'})
def lb2n(l,b):
    l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
def R(t,r): return rotate_xyz(np.eye(3),t,r).T
def theta(n,t,r):
    m=n@R(t,r);psi=np.degrees(np.arccos(np.clip(m[...,0],-1,1)))
    return np.where(m[...,1]>=0,psi,360-psi)
STEPS=[('Untilted',0,0),('+15° tilt',15,0),('+15° tilt + 80.8° roll',15,80.8)]
from matplotlib.colors import ListedColormap
cm=ListedColormap(['#9ecae1','#fcae91'])
ls=np.linspace(-180,180,721);bs=np.linspace(-90,90,361);LL,BB=np.meshgrid(ls,bs);N=lb2n(LL,BB)
fig=plt.figure(figsize=(22,30),facecolor='white');gs=fig.add_gridspec(6,3,width_ratios=[1,2.2,1]);rows=[];score={}
for i,(nm,t,r) in enumerate(STEPS):
  for j,k in enumerate('AB'):
    rr=i*2+j;spots=[];
    for (s,l,b,typ) in SPOTS[k]:
        tt=float(theta(lb2n(l,b),t,r));z=zone(k,tt);m=(z>0)==(typ=='hot');spots.append((s,l,b,typ,tt,z,m))
        rows.append((nm,k,s,typ,l,b,tt,'MAX' if z>0 else 'MIN',m))
    score[(nm,k)]=sum(x[6] for x in spots)
    hdr=f"{nm} — {k} residual zones ({score[(nm,k)]}/{len(spots)})"
    # 1 side circle (model frame: +x anti-CoE right, scan angle θ)
    ax=fig.add_subplot(gs[rr,0]);c,sg=ZS[k]
    for q in range(4):
        lo,hi=c[q],c[q]+((c[(q+1)%4]-c[q])%360);th=np.radians(np.linspace(lo,hi,100))
        ax.fill(np.r_[0,np.cos(th),0],np.r_[0,np.sin(th),0],color=cm(1 if sg[q]>0 else 0),lw=0)
    for cc in c: ax.plot([0,np.cos(np.radians(cc))],[0,np.sin(np.radians(cc))],'k',lw=1.2);ax.text(1.12*np.cos(np.radians(cc)),1.12*np.sin(np.radians(cc)),f"{cc}°",ha='center',va='center',fontsize=8)
    ax.add_patch(plt.Circle((0,0),1,fill=False,lw=1.5))
    for (s,l,b,typ,tt,z,m) in spots:
        a=np.radians(tt);ax.plot(.82*np.cos(a),.82*np.sin(a),'o',mfc='red' if typ=='hot' else 'blue',mec='k',ms=8)
    ax.annotate('anti-CoE (+x)',(1,0),xytext=(1.05,-.18),fontsize=8);ax.text(-1.32,-.18,'CoE (−x)',fontsize=8);ax.plot(0,0,'k^')
    ax.set_xlim(-1.4,1.4);ax.set_ylim(-1.3,1.3);ax.set_aspect('equal');ax.axis('off');ax.set_title("model frame, θ′ from anti-CoE",fontsize=9)
    # 2 Mollweide
    ax=fig.add_subplot(gs[rr,1],projection='mollweide');Z=zone(k,theta(N,t,r))
    ax.pcolormesh(-np.radians(LL),np.radians(BB),Z,cmap=cm,vmin=-1,vmax=1,shading='auto')
    ax.contour(-np.radians(LL),np.radians(BB),Z,[0],colors='k',linewidths=1.2)
    an=R(t,r)[:,0];al=np.degrees(np.arctan2(an[1],an[0]));ab=np.degrees(np.arcsin(an[2]))
    ax.plot(-np.radians(((al+180)%360)-180),np.radians(ab),'k*',ms=12)
    for (s,l,b,typ,tt,z,m) in spots:
        x=-np.radians(((l+180)%360)-180);ax.plot(x,np.radians(b),'o',mfc='red' if typ=='hot' else 'blue',mec='k',ms=8,zorder=5)
        ax.text(x,np.radians(b)+.12,('✓' if m else '✗'),color='green' if m else 'crimson',fontsize=13,weight='bold',ha='center',zorder=6)
    ax.set_xticks(np.radians([-120,-60,0,60,120]));ax.set_xticklabels(['120','60','0','300','240'],fontsize=7);ax.tick_params(labelsize=7);ax.grid(alpha=.3)
    ax.set_title(hdr+"\nGalactic ℓ (increases LEFT)",fontsize=11,weight='bold')
    # 3 globe, orthographic centred on (ℓ=0,b=+20) viewed from outside
    ax=fig.add_subplot(gs[rr,2]);vc=lb2n(-30,20);e1=np.cross([0,0,1],vc);e1/=np.linalg.norm(e1);e2=np.cross(vc,e1)
    vis=(N@vc)>0;X=-(N@e1);Y=N@e2
    ax.scatter(X[vis][::3],Y[vis][::3],c=Z[vis][::3],cmap=cm,vmin=-1,vmax=1,s=2,lw=0)
    ax.add_patch(plt.Circle((0,0),1,fill=False,lw=1.2))
    for (s,l,b,typ,tt,z,m) in spots:
        n=lb2n(l,b)
        if n@vc>0: ax.plot(-(n@e1),n@e2,'o',mfc='red' if typ=='hot' else 'blue',mec='k',ms=8)
    if an@vc>0: ax.plot(-(an@e1),an@e2,'k*',ms=12)
    ax.set_aspect('equal');ax.axis('off');ax.set_title("globe (front hemisphere, ℓ≈330°, b=+20°)",fontsize=9)
lines=[f"{nm} | {k} | {s} ({l},{b}) {typ} → θ′={tt:.0f}° {zz} {'✓' if m else '✗'}" for nm,k,s,typ,l,b,tt,zz,m in rows]
tot=[f"{nm}: A {score[(nm,'A')]}/5, B {score[(nm,'B')]}/7, total {score[(nm,'A')]+score[(nm,'B')]}/12" for nm,_,_ in STEPS]
fig.text(.5,.012,"Bands = A/B residual zones from crossings (A 22/99/262/310, B 14/104/266/321): salmon = MAX (>0), blue = MIN (<0); flat zones, no density/opacity gradient.\n"
 "● red = hot spot, ● blue = cold spot, fixed at published (ℓ,b). ✓ = hot in MAX or cold in MIN.  ★ = anti-CoE axis (physics.rotate_xyz; tilt β, roll γ about the axis).\n"
 "θ′ = arccos(m_x) with m = Rᵀn̂; side-of-ring rule: θ′ kept if model y′ ≥ 0, else 360°−θ′.\n"+"   |   ".join(tot),ha='center',fontsize=11)
plt.subplots_adjust(top=.98,bottom=.05,hspace=.3,wspace=.05)
plt.savefig('rot_steps_nograd.png',dpi=80,facecolor='white')
print("\n".join(lines));print(tot)
