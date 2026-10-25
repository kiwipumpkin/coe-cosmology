import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
sys.path.insert(0,'/workspace/figs_letter_exact_python/makers');sys.path.insert(0,'/tmp')
from physics import rotate_xyz
from spots import SPOTS
M=json.load(open('/tmp/M25.json'));L=np.array(M['l']);RES={k:np.array(M[k]) for k in 'AB'}
def interp(k,t):
    l=np.r_[L-360,L,L+360];v=np.tile(RES[k],3);return np.interp(np.asarray(t)%360,l,v)
def lb2n(l,b):
    l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
def Rmat(t,r): return rotate_xyz(np.eye(3),t,r).T   # columns = model axes in Galactic; R@model=gal
def theta(n,t,r):
    m=n@Rmat(t,r)        # model coords = R^T n
    psi=np.degrees(np.arccos(np.clip(m[...,0],-1,1)))
    return np.where(m[...,1]>=0,psi,360-psi)
STEPS=[('Untilted',0,0),('+15° tilt',15,0),('+15° tilt, +80.8° roll',15,80.8)]
# check rotate_xyz maps +x to anti-CoE (0,+15)
print(Rmat(15,0)[:,0], Rmat(15,80.8)[:,0])
ls=np.linspace(-180,180,721);bs=np.linspace(-90,90,361);LL,BB=np.meshgrid(ls,bs);N=lb2n(LL,BB)
fig,axs=plt.subplots(3,2,figsize=(16,15),subplot_kw={'projection':'mollweide'})
rows=[];score=[]
for i,(nm,t,r) in enumerate(STEPS):
    th=theta(N,t,r);ok=0;tot=0
    for j,k in enumerate('AB'):
        ax=axs[i,j];Z=interp(k,th);vm=np.abs(RES[k]).max()
        ax.pcolormesh(-np.radians(LL),np.radians(BB),Z,cmap='RdBu_r',vmin=-vm,vmax=vm,shading='auto')
        ax.contour(-np.radians(LL),np.radians(BB),Z,[0],colors='k',linewidths=1.5)
        an=Rmat(t,r)[:,0];al=np.degrees(np.arctan2(an[1],an[0]));ab=np.degrees(np.arcsin(an[2]))
        ax.plot(-np.radians(((al+180)%360)-180),np.radians(ab),'k*',ms=14);ax.plot(-np.radians(((al)%360)-180),np.radians(-ab),'kX',ms=11)
        for (s,l,b,typ) in SPOTS[k]:
            tt=float(theta(lb2n(l,b),t,r));v=float(interp(k,tt));m=(v>0)==(typ=='hot');ok+=m;tot+=1
            rows.append((nm,k,s,typ,l,b,tt,v,m))
            x=-np.radians(((l+180)%360)-180)
            ax.plot(x,np.radians(b),'o',mfc='red' if typ=='hot' else 'blue',mec='k' if m else 'yellow',mew=2,ms=10,zorder=5,clip_on=False)
            ax.annotate(f"{s}\nθ'={tt:.0f}° {v:+.1f}% {'✓' if m else '✗'}",(x,np.radians(b)),xytext=(4,5),textcoords='offset points',fontsize=7.5,bbox=dict(fc='white',alpha=.7,lw=0))
        ax.set_xticks(np.radians([-150,-120,-90,-60,-30,0,30,60,90,120,150]));ax.set_xticklabels(['150','120','90','60','30','0','330','300','270','240','210'],fontsize=7)
        ax.grid(alpha=.3)
        sc=sum(1 for rr in rows if rr[0]==nm and rr[1]==k and rr[8])
        ax.set_title(f"{nm} — residual {k} (score {sc}/{len(SPOTS[k])})",fontsize=12,pad=12)
    score.append((nm,ok,tot))
fig.text(.5,.035,"Colour: residual % (red >0, blue <0); black line = zero contour. ● red = hot spot, ● blue = cold spot; black edge = match, yellow edge = miss.\n★ anti-CoE (+x model axis), ✖ CoE. θ' = scan angle from anti-CoE, cosψ'=[R_roll(γ)R_tilt(β)D̂]·n̂, side from model y' (physics.rotate_xyz). Galactic ℓ increases left.",ha='center',fontsize=10)
fig.text(.5,.01,"   ".join(f"{n}: {o}/{t}" for n,o,t in score),ha='center',fontsize=13,weight='bold')
plt.subplots_adjust(bottom=.08,top=.96,hspace=.25,wspace=.08)
plt.savefig('rot_steps_deckgeom.png',dpi=110,facecolor='white')
for rr in rows: print(f"{rr[0]:24s} {rr[1]} {rr[2]:18s} {rr[3]:4s} ({rr[4]},{rr[5]}) θ'={rr[6]:6.1f} {rr[7]:+5.2f} {'✓' if rr[8] else '✗'}")
print(score)
