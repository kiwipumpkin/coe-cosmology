import sys,json,numpy as np,matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
sys.path.insert(0,'/workspace/figs_letter_exact_python/makers'); sys.path.insert(0,'/tmp')
from physics import rotate_xyz,TILT_DEG,ROLL_DEG
from spots import SPOTS
d=json.load(open('/tmp/M25.json')); L=np.array(d['l'])
R=rotate_xyz(np.eye(3),tilt_deg=TILT_DEG,roll_deg=ROLL_DEG)  # v_gal=v_mod@R
def field(key,l,b):
    l=np.radians(l);b=np.radians(b)
    v=np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
    m=v@R.T; lm=np.degrees(np.arctan2(m[...,1],m[...,0]))%360
    return np.interp(lm,L,np.array(d[key]),period=360)
def oldzone(k,l):
    hiA=(22<=l<=99) or (262<=l<=310); loB=(14<=l<=104) or (266<=l<=321)
    return (1 if hiA else -1) if k=='A' else (-1 if loB else 1)
lg=np.linspace(-180,180,721); bg=np.linspace(-90,90,361); LG,BG=np.meshgrid(lg,bg)
fig=plt.figure(figsize=(16,24),facecolor='w')
rows=[];sc={'A':[0,0],'B':[0,0]}
for i,k in enumerate('AB'):
    ax=fig.add_axes([0.06,0.67-0.33*i,0.88,0.28],projection='mollweide')
    Z=field(k,LG%360,BG); lim=np.abs(Z).max()
    cf=ax.pcolormesh(-np.radians(LG),np.radians(BG),Z,cmap='RdBu_r',norm=TwoSlopeNorm(0,-lim,lim),shading='auto')
    ax.contour(-np.radians(LG),np.radians(BG),Z,[0],colors='k',linewidths=2.5)
    ax.set_xticklabels([f'{int(t)}°' for t in [150,120,90,60,30,0,330,300,270,240,210]],fontsize=12)
    ax.grid(color='0.5',lw=0.5)
    ax.set_title(f'{k} residual (%), transformed: +{TILT_DEG:g}° tilt, +{ROLL_DEG:g}° roll; spots at published (ℓ,b), unrotated',fontsize=16)
    for n,l,b,hc in SPOTS[k]:
        v=float(field(k,l,b)); ok=(v>0)==(hc=='hot'); oo=(oldzone(k,l)>0)==(hc=='hot')
        sc[k][0]+=ok;sc[k][1]+=1; rows.append([n,f'{l:g}',f'{b:+g}',k,hc,f'{v:+.1f}%','✓' if ok else '✗','✓' if oo else '✗'])
        x=-np.radians(((l+180)%360)-180)
        ax.plot(x,np.radians(b),'o',ms=18,mfc='red' if hc=='hot' else 'blue',mec='k',mew=2)
        ax.annotate(n,(x,np.radians(b)),xytext=(0,16),textcoords='offset points',ha='center',fontsize=13,weight='bold',color='red' if hc=='hot' else 'blue',bbox=dict(fc='w',ec='none',alpha=.8))
    cb=fig.colorbar(cf,ax=ax,orientation='horizontal',fraction=0.05,pad=0.06); cb.set_label(f'{k}: 100·(data−model)/data, ±25° mean; black = zero contour',fontsize=13)
ax=fig.add_axes([0.04,0.02,0.92,0.25]);ax.axis('off')
t=ax.table(cellText=rows,colLabels=['Spot','ℓ','b','Tied to','Hot/cold','Transformed map','2D match','Old ℓ-only'],loc='upper center',cellLoc='center')
t.auto_set_font_size(False);t.set_fontsize(14);t.scale(1,2.3);t.auto_set_column_width(range(8))
tot=sc['A'][0]+sc['B'][0]
ax.text(0.5,0.0,f"2D score: A {sc['A'][0]}/{sc['A'][1]}, B {sc['B'][0]}/{sc['B'][1]}, total {tot}/12   (old ℓ-only: 11/12)",ha='center',fontsize=18,weight='bold',transform=ax.transAxes)
fig.savefig('combined_rot2d.png',dpi=110,facecolor='w')
for r in rows:print(r)
print(sc,R)
