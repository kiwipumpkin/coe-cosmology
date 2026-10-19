import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
G=42.0;Gy=G*3.26156;R=0.5;D={k:np.load(f'z4_{k}.npz') for k in ('nomass','m02')}
lo=min(D[k]['o'][:,2].min() for k in D)*Gy;hi=max(D[k]['o'][:,2].max() for k in D)*Gy;norm=plt.Normalize(lo,hi);cm='plasma'
fig,axs=plt.subplots(2,2,figsize=(20,16),gridspec_kw=dict(width_ratios=[1,1.1]));th=np.linspace(0,2*np.pi,721)
for r,(k,lab) in enumerate((('nomass','no mass'),('m02','0.2 of critical'))):
    d=D[k];o=d['o'];re,rt,phi=d['re'],d['rt'],d['phi'];xo,xc,t0,tc=[float(d[q]) for q in('xo','xc','t0','tc')];rv=t0-tc;te=o[:,2]*Gy
    print(f'{lab}: emission time {te.min():.3f}-{te.max():.3f} Gyr (0° {te[0]:.3f}, 90° {te[90]:.3f}, 180° {te[180]:.3f})')
    P=np.r_[phi,-phi[::-1]];TE=np.r_[te,te[::-1]];RE=np.r_[re,re[::-1]]*G;RT=np.r_[rt,rt[::-1]]*G;X0=xo*G;dt=np.hypot(rt*np.cos(phi)-xo,rt*np.sin(phi))
    for c_,ax in enumerate(axs[r]):
        ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
        for RR,lw in ((RT,4),(RE,5)):
            x,y=RR*np.cos(P),RR*np.sin(P);seg=np.stack([np.c_[x[:-1],y[:-1]],np.c_[x[1:],y[1:]]],1)
            lc=LineCollection(seg,cmap=cm,norm=norm,lw=lw);lc.set_array(0.5*(TE[:-1]+TE[1:]));ax.add_collection(lc)
        ax.plot((xo+rv*np.cos(th))*G,rv*np.sin(th)*G,'g-',lw=2)
        for s in (1,-1):
            for i in range(0,181,15):
                cc,sn=np.cos(s*phi[i]),np.sin(s*phi[i]);x1,y1,x2,y2=re[i]*G*cc,re[i]*G*sn,rt[i]*G*cc,rt[i]*G*sn
                ax.plot([0,x1],[0,y1],'--',color='gray',lw=.7);ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='gray',lw=.9));ax.plot([x2,X0],[y2,0],color='darkorange',lw=1.1)
        ax.plot(X0,0,'k*',ms=15);ax.annotate('us today',(X0,0),(6,8),textcoords='offset points')
        ax.plot(xc*G,0,'ko',mfc='none',ms=8);ax.annotate('us at crossing',(xc*G,0),(6,-14),textcoords='offset points')
        ax.plot(0,0,'k+',ms=18,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points');ax.set_aspect('equal');ax.grid(alpha=.3);ax.set_xlabel('Gpc (anti-CoE side to the right)');ax.set_ylabel('Gpc')
        if c_==0: ax.set_xlim(-R*G*1.1,RT.max()*np.cos(0)+3);ax.set_ylim(-R*G*1.1,R*G*1.1);ax.set_title(f'{lab}, full view; our speed 4 c; today {t0*Gy:.1f} Gyr; green {rv*G:.2f} Gpc',fontsize=12)
        else:
            X=np.r_[RE*np.cos(P),RT*np.cos(P)];Y=np.r_[RE*np.sin(P),RT*np.sin(P)];ax.set_xlim(X.min()-1,X.max()+1);ax.set_ylim(Y.min()-1,Y.max()+1)
            ax.set_title(f'{lab}, zoom; straight red-to-us: 0° {dt[0]*G:.2f}, 90° {dt[90]*G:.2f}, 180° {dt[180]*G:.2f} Gpc',fontsize=12)
h=[plt.Line2D([],[],color=plt.get_cmap(cm)(0.6),lw=4),plt.Line2D([],[],color='g',lw=2),plt.Line2D([],[],color='darkorange'),plt.Line2D([],[],color='gray',ls='--')]
fig.legend(h,['thick coloured line far out = same emitters today (was red); thick coloured arc on sphere = emission points; colour = recombination time','green: c × time since we crossed, centred on us today','orange: straight line from today position to us (every 15°)','gray: dashed radial line from CoE, arrow = displacement emission → today'],loc='lower left',bbox_to_anchor=(0.05,0.92),fontsize=12,frameon=False)
sm=plt.cm.ScalarMappable(norm=norm,cmap=cm);cb=fig.colorbar(sm,ax=axs,shrink=0.6,pad=0.02);cb.set_label('recombination time (Gyr)',fontsize=13)
fig.savefig('emit_cross_4c_2x2_coloured.png',dpi=90,bbox_inches='tight')
