# rev7 regeneration of rev6 Figs 26-29 (sky-residual chain), one script, everything computed:
# (1) recon_n3.py prefix (cosmology/figures/fourier_resid): Migkas A,B data, locked CoE fits (l50/curves.json), shared n=1-3 harmonics removed;
# (2) 25-deg running mean of the % residual (was /tmp/M25.json via dumpM.py); (3) zero crossings computed from (2) (rev6 used hand-entered
# /tmp/R.json windows); (4) lift 15 deg + roll 45 deg about the CoE axis (sky_mapA_al.py); (5) residual at each published spot -> table and score.
# The locked Migkas fit is unchanged in rev7 (D=34.96 refit has the same extrema 20/200/127/272), so these are the rev7 values.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import sys;sys.path.insert(0,'/workspace/cosmology/figures');H=os.path.dirname(os.path.abspath(__file__));os.chdir('/workspace/cosmology/figures/fourier_resid')
exec(open('recon_n3.py').read().split("cols=['purple'")[0]);exec(open('/workspace/coepdf/book/recovered/spots.py').read())
res={};MM={'l':list(map(float,l))}
for nm,Y,mb,k in (('A',A0,AFl,0),('B',B0,BLl,1)):
    m=mb.copy()
    for n,p in steps2[:3]:m=m+T(n,p,th,k)
    res[nm]=100*(Y-m)/Y;MM[nm]=[float(np.mean(res[nm][np.abs((l-c+180)%360-180)<=25])) for c in l]
gq=np.arange(0,360,0.25);Z={}
for nm in 'AB':
    G=np.interp(gq,MM['l'],MM[nm],period=360);Z[nm]=[float(gq[i]) for i in range(len(gq)) if np.sign(G[i])!=np.sign(G[(i+1)%len(gq)])]
json.dump(dict(M25=MM,zero_crossings=Z),open(os.path.join(H,'R26_29_sky_chain.json'),'w'))
TILT,ROLL=-15.0,45.0
def vec(l_,b_):l_,b_=np.radians(l_),np.radians(b_);return np.stack([np.cos(b_)*np.cos(l_),np.cos(b_)*np.sin(l_),np.sin(b_)],-1)
def rot(ax,a):
    ax=ax/np.linalg.norm(ax);a=np.radians(a);K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]]);return np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K
R1=rot(np.array([0,1.,0]),TILT);coe=R1@vec(180,0);RT=rot(coe,ROLL)@R1
def lm(l_,b_):m=vec(l_,b_)@RT;return np.degrees(np.arctan2(m[...,1],m[...,0]))%360
def moll(l_,b_):
    lam=np.radians(((np.asarray(l_,float)+180)%360)-180);b_=np.radians(np.asarray(b_,float));t=b_.copy()
    for _ in range(60):t=t-(2*t+np.sin(2*t)-np.pi*np.sin(b_))/(2+2*np.cos(2*t)+1e-12)
    return -2*np.sqrt(2)/np.pi*lam*np.cos(t),np.sqrt(2)*np.sin(t)
MK=['o','s','D','^','v','p','h'];rows=[];score={}
# Fig 26: residual vs model-frame angle
fig,axs=plt.subplots(2,1,figsize=(18,11))
def rmn(r,w):return np.array([np.mean(r[np.abs((l-c+180)%360-180)<=w]) for c in l])
for a,nm in zip(axs,'AB'):
    a.axhline(0,color='grey');h=[]
    for w,c in zip((15,25,35),plt.cm.viridis(np.linspace(0,.7,3))):
        M=rmn(res[nm],w);i,j=np.argmax(M),np.argmin(M);h.append(a.plot(l,M,color=c,lw=2.2,label=f'±{w}° (max {M[i]:+.1f}% at {l[i]:.0f}°, min {M[j]:+.1f}% at {l[j]:.0f}°)')[0])
    yt=a.get_ylim()[1]
    for i,(t,sl,sb,c) in enumerate(SPOTS[nm]):
        x=float(lm(sl,sb));a.axvline(x,color=CC[c],ls='-.');a.plot(x,yt*.9,MK[i],ms=14,mfc=CC[c],mec='k');h.append(Line2D([],[],ls='',marker=MK[i],mfc=CC[c],mec='k',ms=10,label=f'{i+1} {t}'))
    for z in Z[nm]:a.axvline(z,color='k',ls='--')
    a.set_xlim(0,360);a.set_xlabel('model-frame angle [deg]');a.set_ylabel(f'{nm} running-mean residual [%]');a.legend(handles=h,ncol=3,fontsize=9,loc='upper center',bbox_to_anchor=(.5,-.12),frameon=False)
fig.tight_layout();fig.savefig(os.path.join(H,'R26_sky_angle.png'),dpi=90)
# Figs 27/28: residual maps
X,Y_=np.meshgrid(np.linspace(-2*np.sqrt(2),2*np.sqrt(2),1000),np.linspace(-np.sqrt(2),np.sqrt(2),500));ok=(X/(2*np.sqrt(2)))**2+(Y_/np.sqrt(2))**2<=1
t_=np.arcsin(np.clip(Y_/np.sqrt(2),-1,1));Bg=np.degrees(np.arcsin(np.clip((2*t_+np.sin(2*t_))/np.pi,-1,1)));Lg=np.degrees(-np.pi*X/(2*np.sqrt(2)*np.cos(t_)+1e-12))%360;LM=lm(Lg,Bg)
for nm,fn in (('A','R27_mapA.png'),('B','R28_mapB.png')):
    fig,ax=plt.subplots(figsize=(15,8.5));v=np.interp(LM,MM['l'],MM[nm],period=360);gm=np.abs(v[ok]).max()
    ax.imshow(np.where(ok,v,np.nan),extent=(X.min(),X.max(),Y_.min(),Y_.max()),origin='lower',cmap='bwr',vmin=-gm,vmax=gm,alpha=.7)
    for z in Z[nm]:ax.contour(X,Y_,np.where(ok,(LM-z+180)%360-180,np.nan),[0],colors='k',linewidths=2)
    tt=np.linspace(0,2*np.pi,400);ax.plot(2*np.sqrt(2)*np.cos(tt),np.sqrt(2)*np.sin(tt),'k');sc=0;h=[]
    for i,(t,sl,sb,c) in enumerate(SPOTS[nm]):
        val=float(np.interp(lm(sl,sb),MM['l'],MM[nm],period=360));okk=(val>0)==(c=='hot');amb=abs(val)<0.1;sc+=okk and not amb
        rows.append((nm,t,f'{sl:g}, {sb:+g}',c,'> 0' if c=='hot' else '< 0',f'{val:+.2f}','?' if amb else ('✓' if okk else '✗')))
        xs,ys=map(float,moll(sl,sb));ax.plot(xs,ys,MK[i],ms=24,mfc=CC[c],mec='k',mew=2);ax.text(xs,ys,str(i+1),ha='center',va='center',color='w',fontweight='bold')
        ax.annotate(f'{t} {val:+.2f}%',(xs,ys),xytext=(12,12),textcoords='offset points',fontsize=11,color=CC[c],bbox=dict(fc='w',ec=CC[c]))
    score[nm]=(sc,len(SPOTS[nm]));ax.set_title(f'{nm} residual (lift 15°, roll 45° about CoE axis): {sc}/{len(SPOTS[nm])} match',fontsize=14);ax.axis('off');ax.set_aspect('equal')
    fig.tight_layout();fig.savefig(os.path.join(H,fn),dpi=90)
# Fig 29: table, values computed above
fig,ax=plt.subplots(figsize=(15,5.2));ax.axis('off');tb=ax.table(cellText=rows,colLabels=['Tied to','Spot','l, b','Hot/cold','Needs','Residual (%)','Match'],loc='center',cellLoc='center')
tb.auto_set_font_size(False);tb.set_fontsize(11);tb.scale(1,1.5)
for i,r in enumerate(rows,1):tb[i,6].get_text().set_color({'✓':'green','?':'darkorange'}.get(r[6],'red'))
ax.set_title(f'Score (computed): A {score["A"][0]}/{score["A"][1]}, B {score["B"][0]}/{score["B"][1]}',fontsize=14);fig.savefig(os.path.join(H,'R29_spot_tab.png'),dpi=110,bbox_inches='tight')
json.dump(dict(M25=MM,zero_crossings=Z,rows=rows,score=score),open(os.path.join(H,'R26_29_sky_chain.json'),'w'),indent=1);print(Z,score);[print(r) for r in rows]
