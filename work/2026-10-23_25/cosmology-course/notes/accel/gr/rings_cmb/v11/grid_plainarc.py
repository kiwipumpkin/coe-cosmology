import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
cfg=json.loads(sys.argv[1]);G=cfg['D'];Gy=G*3.26156;R=0.5;rows=[]
for r in cfg['rows']:
    if r['type']=='traced':
        d=np.load(r['file']);o=d['o'];ok=np.isfinite(o[:,1]);GG=(4.2/(float(d['t0'])-float(d['tc'])) if r.get('green42') else cfg['D'])
        rows.append(dict(D=(4.2/(float(d['t0'])-float(d['tc'])) if r.get('green42') else cfg['D']),lab=r['lab'],re=d['re'],rt=d['rt'],phi=d['phi'],te=o[:,2]*GG*3.26156,ok=ok,xo=float(d['xo']),xc=float(d['xc']),t0=float(d['t0']),tc=float(d['tc'])))
        if (~ok).any(): print(r['lab'],'UNRESOLVED',np.arange(181)[~ok].tolist())
    else:
        u=r['u'];t0=1/u;tc=R/u;a=np.radians(np.arange(0,181,1.0));rho=0.4*R;xc_,yc_=u*tc+rho*np.cos(a),rho*np.sin(a)
        phi=np.arctan2(yc_,xc_);rc=np.hypot(xc_,yc_)
        rows.append(dict(lab=r['lab'],re=np.full(181,R),rt=rc*t0/tc,phi=phi,te=np.full(181,tc*Gy),ok=np.ones(181,bool),xo=u*t0,xc=u*tc,t0=t0,tc=tc,sim=dict(rc=rc)))
lo=min(x['te'][x['ok']].min() for x in rows);hi=max(x['te'][x['ok']].max() for x in rows);norm=plt.Normalize(lo,hi);cm='plasma'
n=len(rows);AR=1.5;pw,ph=7.5,7.5/AR;mx,my,top=0.9,1.0,2.1
fig=plt.figure(figsize=(2*pw+3*mx+1.6,n*ph+(n)*my+top))
FW,FH=fig.get_size_inches();axs=np.array([[fig.add_axes([(mx+c*(pw+mx))/FW,1-(top+(r+1)*ph+r*my)/FH,pw/FW,ph/FH]) for c in range(2)] for r in range(n)])
th=np.linspace(0,2*np.pi,721)
cax=fig.add_axes([(2*pw+3*mx-0.3)/FW,0.15,0.25/FW,0.6])
def fit(ax,x0,x1,y0,y1):
    cx,cy,w,h=(x0+x1)/2,(y0+y1)/2,x1-x0,y1-y0
    if w/h<AR: w=h*AR
    else: h=w/AR
    ax.set_xlim(cx-w/2,cx+w/2);ax.set_ylim(cy-h/2,cy+h/2);th=np.linspace(0,2*np.pi,721)
for k,x in enumerate(rows):
    G=x.get('D',cfg['D']);Gy=G*3.26156
    re,rt,phi,te,ok=x['re'][x['ok']],x['rt'][x['ok']],x['phi'][x['ok']],x['te'][x['ok']],x['ok'];xo,xc,t0,tc=x['xo'],x['xc'],x['t0'],x['tc'];rv=t0-tc;X0=xo*G
    dt=np.hypot(x['rt']*np.cos(x['phi'])-xo,x['rt']*np.sin(x['phi']))
    print(f"{x['lab']}: recomb time {te.min():.3f}-{te.max():.3f} Gyr; today {t0*Gy:.3f} Gyr; our crossing {tc*Gy:.3f}; green {rv*G:.3f} Gpc; red-to-us 0/90/180 {dt[0]*G:.3f}/{dt[90]*G:.3f}/{dt[180]*G:.3f} Gpc; arc half-width {np.degrees(np.abs(phi).max()):.3f} deg; green inside red {np.nanmin(dt)>=rv}")
    P=np.r_[phi,-phi[::-1]];TE=np.r_[te,te[::-1]];RE=np.r_[re,re[::-1]]*G;RT=np.r_[rt,rt[::-1]]*G
    for c_,ax in enumerate(axs[k]):
        ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
        xx,yy=RT*np.cos(P),RT*np.sin(P);seg=np.stack([np.c_[xx[:-1],yy[:-1]],np.c_[xx[1:],yy[1:]]],1)
        lc=LineCollection(seg,cmap=cm,norm=norm,lw=5,zorder=3);lc.set_array(0.5*(TE[:-1]+TE[1:]));ax.add_collection(lc)
        ax.plot(RE*np.cos(P),RE*np.sin(P),'-',color='blue',lw=5,solid_capstyle='butt',zorder=3)
        ax.plot((xo+rv*np.cos(th))*G,rv*np.sin(th)*G,'g-',lw=2)
        for s in (1,-1):
            for i in range(0,181,15):
                if not ok[i]: continue
                cc,sn=np.cos(s*x['phi'][i]),np.sin(s*x['phi'][i]);x1,y1,x2,y2=x['re'][i]*G*cc,x['re'][i]*G*sn,x['rt'][i]*G*cc,x['rt'][i]*G*sn
                ax.plot([0,x1],[0,y1],'--',color='gray',lw=.7);ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='gray',lw=.9));ax.plot([x2,X0],[y2,0],color='darkorange',lw=1.1)
        ax.plot(X0,0,'k*',ms=14);ax.annotate('us today',(X0,0),(6,8),textcoords='offset points')
        ax.plot(xc*G,0,'ko',mfc='none',ms=8);ax.annotate('us at crossing',(xc*G,0),(6,-14),textcoords='offset points')
        ax.plot(0,0,'k+',ms=16,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points');ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.3);ax.set_ylabel('Gpc')
        if k==n-1: ax.set_xlabel('Gpc (anti-CoE side to the right)')
        if c_==0:
            yl=max(R*G*1.1,np.abs(RT*np.sin(P)).max()*1.1);fit(ax,-R*G*1.1,RT.max()+0.05*G,-yl,yl)
            ax.set_title(f"{x['lab']}\nfull view; today {t0*Gy:.1f} Gyr; green {rv*G:.2f} Gpc",fontsize=11,loc='left')
        else:
            XX=np.r_[RE*np.cos(P),RT*np.cos(P)];YY=np.r_[RE*np.sin(P),RT*np.sin(P)];m=0.03*(XX.max()-XX.min());fit(ax,XX.min()-m,XX.max()+m,YY.min()-m,YY.max()+m)
            ax.set_title(f"zoom; today position to us:\n0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc",fontsize=11,loc='left')
h=[plt.Line2D([],[],color='blue',lw=5),plt.Line2D([],[],color=plt.get_cmap(cm)(0.6),lw=5),plt.Line2D([],[],color='g',lw=2),plt.Line2D([],[],color='darkorange'),plt.Line2D([],[],color='gray',ls='--')]
fig.legend(h,['blue arc on recombination sphere: emission points','big bubble: same emitters today, coloured by recombination time','green: c × time since we crossed, centred on us today','orange: straight line from today position to us (every 15°)','gray: dashed radial line from CoE, arrow = displacement emission → today'],loc='lower left',bbox_to_anchor=(0.02,1-1.5/FH),fontsize=12,frameon=False,ncol=2)
fig.suptitle(cfg['title'],y=1-0.25/FH,fontsize=13)
sm=plt.cm.ScalarMappable(norm=norm,cmap=cm);cb=fig.colorbar(sm,cax=cax);cb.set_label('recombination time (Gyr)',fontsize=13)
fig.savefig(cfg['out'],dpi=90)
print('axes sizes',set(tuple(np.round(a.get_position().bounds[2:],4)) for a in axs.ravel()))
