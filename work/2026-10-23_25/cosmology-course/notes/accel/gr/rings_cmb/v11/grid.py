import sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
cfg=json.loads(sys.argv[1]);G=cfg['D'];Gy=G*3.26156;R=0.5;rows=[]
for r in cfg['rows']:
    if r['type']=='traced':
        d=np.load(r['file']);o=d['o'];ok=np.isfinite(o[:,1])
        rows.append(dict(lab=r['lab'],re=d['re'],rt=d['rt'],phi=d['phi'],te=o[:,2]*Gy,ok=ok,xo=float(d['xo']),xc=float(d['xc']),t0=float(d['t0']),tc=float(d['tc'])))
        if (~ok).any(): print(r['lab'],'UNRESOLVED',np.arange(181)[~ok].tolist())
    else:
        u=r['u'];t0=1/u;tc=R/u;a=np.radians(np.arange(0,181,1.0));rho=0.4*R;xc_,yc_=u*tc+rho*np.cos(a),rho*np.sin(a)
        phi=np.arctan2(yc_,xc_);rc=np.hypot(xc_,yc_)
        rows.append(dict(lab=r['lab'],re=np.full(181,R),rt=rc*t0/tc,phi=phi,te=np.full(181,tc*Gy),ok=np.ones(181,bool),xo=u*t0,xc=u*tc,t0=t0,tc=tc,sim=dict(rc=rc)))
lo=min(x['te'][x['ok']].min() for x in rows);hi=max(x['te'][x['ok']].max() for x in rows);norm=plt.Normalize(lo,hi);cm='plasma'
n=len(rows);fig,axs=plt.subplots(n,2,figsize=(19,4.9*n),gridspec_kw=dict(width_ratios=[1,1.15]),squeeze=False);th=np.linspace(0,2*np.pi,721)
for k,x in enumerate(rows):
    re,rt,phi,te,ok=x['re'][x['ok']],x['rt'][x['ok']],x['phi'][x['ok']],x['te'][x['ok']],x['ok'];xo,xc,t0,tc=x['xo'],x['xc'],x['t0'],x['tc'];rv=t0-tc;X0=xo*G
    dt=np.hypot(x['rt']*np.cos(x['phi'])-xo,x['rt']*np.sin(x['phi']))
    print(f"{x['lab']}: recomb time {te.min():.3f}-{te.max():.3f} Gyr; today {t0*Gy:.3f} Gyr; our crossing {tc*Gy:.3f}; green {rv*G:.3f} Gpc; red-to-us 0/90/180 {dt[0]*G:.3f}/{dt[90]*G:.3f}/{dt[180]*G:.3f} Gpc; arc half-width {np.degrees(np.abs(phi).max()):.3f} deg; green inside red {np.nanmin(dt)>=rv}")
    P=np.r_[phi,-phi[::-1]];TE=np.r_[te,te[::-1]];RE=np.r_[re,re[::-1]]*G;RT=np.r_[rt,rt[::-1]]*G
    for c_,ax in enumerate(axs[k]):
        ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
        for RR,lw in ((RT,4),(RE,5)):
            xx,yy=RR*np.cos(P),RR*np.sin(P);seg=np.stack([np.c_[xx[:-1],yy[:-1]],np.c_[xx[1:],yy[1:]]],1)
            lc=LineCollection(seg,cmap=cm,norm=norm,lw=lw);lc.set_array(0.5*(TE[:-1]+TE[1:]));ax.add_collection(lc)
        ax.plot((xo+rv*np.cos(th))*G,rv*np.sin(th)*G,'g-',lw=2)
        for s in (1,-1):
            for i in range(0,181,15):
                if not ok[i]: continue
                cc,sn=np.cos(s*x['phi'][i]),np.sin(s*x['phi'][i]);x1,y1,x2,y2=x['re'][i]*G*cc,x['re'][i]*G*sn,x['rt'][i]*G*cc,x['rt'][i]*G*sn
                ax.plot([0,x1],[0,y1],'--',color='gray',lw=.7);ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='gray',lw=.9));ax.plot([x2,X0],[y2,0],color='darkorange',lw=1.1)
        ax.plot(X0,0,'k*',ms=14);ax.annotate('us today',(X0,0),(6,8),textcoords='offset points')
        ax.plot(xc*G,0,'ko',mfc='none',ms=8);ax.annotate('us at crossing',(xc*G,0),(6,-14),textcoords='offset points')
        ax.plot(0,0,'k+',ms=16,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points');ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_ylabel('Gpc')
        if k==n-1: ax.set_xlabel('Gpc (anti-CoE side to the right)')
        if c_==0:
            ax.set_xlim(-R*G*1.1,RT.max()+0.05*G);yl=max(R*G*1.1,np.abs(RT*np.sin(P)).max()*1.1);ax.set_ylim(-yl,yl)
            ax.set_title(f"{x['lab']}\nfull view; today {t0*Gy:.1f} Gyr; green {rv*G:.2f} Gpc",fontsize=11,loc='left')
        else:
            XX=np.r_[RE*np.cos(P),RT*np.cos(P)];YY=np.r_[RE*np.sin(P),RT*np.sin(P)];m=0.03*(XX.max()-XX.min());ax.set_xlim(XX.min()-m,XX.max()+m);ax.set_ylim(YY.min()-m,YY.max()+m)
            ax.set_title(f"zoom; today position to us:\n0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc",fontsize=11,loc='left')
h=[plt.Line2D([],[],color=plt.get_cmap(cm)(0.6),lw=4),plt.Line2D([],[],color='g',lw=2),plt.Line2D([],[],color='darkorange'),plt.Line2D([],[],color='gray',ls='--')]
fig.legend(h,['coloured: outer line = same emitters today, arc on sphere = emission points; colour = recombination time','green: c × time since we crossed, centred on us today','orange: straight line from today position to us (every 15°)','gray: dashed radial line from CoE, arrow = displacement emission → today'],loc='lower left',bbox_to_anchor=(0.02,0.905),fontsize=12,frameon=False,ncol=2)
fig.suptitle(cfg['title'],y=0.975,fontsize=13)
fig.subplots_adjust(hspace=0.32,wspace=0.12,top=0.86)
sm=plt.cm.ScalarMappable(norm=norm,cmap=cm);cb=fig.colorbar(sm,ax=axs,shrink=0.5,pad=0.015);cb.set_label('recombination time (Gyr)',fontsize=13)
fig.savefig(cfg['out'],dpi=90,bbox_inches='tight')
