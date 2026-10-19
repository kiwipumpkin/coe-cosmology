import pickle,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
B=pickle.load(open('best.pkl','rb'));res=B['res'];R=0.5
AR=1.5;pw,ph=7.5,7.5/AR;mx,my,top=0.9,1.0,2.1;n=2
fig=plt.figure(figsize=(2*pw+3*mx,n*ph+n*my+top));FW,FH=fig.get_size_inches()
axs=np.array([[fig.add_axes([(mx+c*(pw+mx))/FW,1-(top+(r+1)*ph+r*my)/FH,pw/FW,ph/FH]) for c in range(2)] for r in range(n)])
th=np.linspace(0,2*np.pi,721)
def fit(ax,x0,x1,y0,y1):
    cx,cy,w,h=(x0+x1)/2,(y0+y1)/2,x1-x0,y1-y0
    if w/h<AR: w=h*AR
    else: h=w/AR
    ax.set_xlim(cx-w/2,cx+w/2);ax.set_ylim(cy-h/2,cy+h/2)
for k,(tag,lab) in enumerate((('nomass','no mass'),('mass','fitted core mass 0.0154 of critical'))):
    m=res[tag]['m'];r=[x for x in res[tag]['r'] if 'Pe' in x];G=4.2/(m['t0']-m['tc']);Gy=G*3.26156
    Pe=np.array([x['Pe'] for x in r]);rt=np.array([x['rt'] for x in r]);ok=np.array([np.isfinite(x['opz']) and x['opz']>0 for x in r])
    ph_=np.arctan2(Pe[:,1],Pe[:,0]);PA=np.r_[ph_,-ph_[::-1]];RT=np.r_[rt,rt[::-1]]*G
    rv=(m['t0']-m['tc'])*G;X0=G
    print(tag,'D',G,'today',m['t0']*Gy,'our crossing',m['tc']*Gy,'arc half-width deg',np.degrees(np.abs(ph_).max()),'red radius range Gpc',rt.min()*G,rt.max()*G,'light travel Gpc',min(x['s'] for x in r)*G,max(x['s'] for x in r)*G)
    for c_,ax in enumerate(axs[k]):
        ax.fill(R*G*np.cos(th),R*G*np.sin(th),color='lightblue',alpha=.3);ax.plot(R*G*np.cos(th),R*G*np.sin(th),'b-',lw=1)
        ax.plot(RT*np.cos(PA),RT*np.sin(PA),'r-',lw=3,zorder=3)
        ax.plot(R*G*np.cos(PA),R*G*np.sin(PA),'-',color='blue',lw=5,solid_capstyle='butt',zorder=4)
        ax.plot((1+(m['t0']-m['tc'])*np.cos(th))*G,(m['t0']-m['tc'])*np.sin(th)*G,'g-',lw=2)
        idx=[i for i,x in enumerate(r) if abs(round(x['psi'],2)*100%250)<1e-6]
        for s in (1,-1):
            for i in idx:
                cc,sn=np.cos(s*ph_[i]),np.sin(s*ph_[i]);x1,y1=R*G*cc,R*G*sn;x2,y2=rt[i]*G*cc,rt[i]*G*sn
                ax.plot([0,x1],[0,y1],'--',color='gray',lw=.7);ax.annotate('',(x2,y2),(x1,y1),arrowprops=dict(arrowstyle='->',color='gray',lw=.8))
                ax.plot([x1,X0],[y1,0],color='darkorange',lw=1.2,ls='-' if ok[i] else ':')
        bad=~ok
        if bad.any():
            for s in (1,-1): ax.plot(R*G*np.cos(s*ph_[bad]),R*G*np.sin(s*ph_[bad]),'m.',ms=5,zorder=5)
        ax.plot(X0,0,'k*',ms=14,zorder=6);ax.annotate('us today',(X0,0),(6,8),textcoords='offset points')
        ax.plot(m['tc']*0+R*G,0,'ko',mfc='none',ms=8,zorder=6);ax.annotate('us at crossing',(R*G,0),(6,-14),textcoords='offset points')
        ax.plot(0,0,'k+',ms=16,mew=2);ax.annotate('CoE',(0,0),(-20,10),textcoords='offset points');ax.set_aspect('equal',adjustable='box');ax.grid(alpha=.3);ax.set_ylabel('Gpc')
        if k==n-1: ax.set_xlabel('Gpc (anti-CoE side to the right)')
        if c_==0:
            yl=max(R*G*1.1,np.abs(RT*np.sin(PA)).max()*1.1);fit(ax,-R*G*1.1,RT.max()*1.03,-yl,yl)
            ax.set_title(f"Static space, our speed 1.052 c, {lab}\nfull view; today {m['t0']*Gy:.2f} Gyr; green {rv:.2f} Gpc; D {G:.2f} Gpc",fontsize=11,loc='left')
        else:
            fit(ax,-0.05*G,(1+(m['t0']-m['tc']))*G*1.05,-0.62*G,0.62*G)
            ax.set_title(f"zoom; emission point to us: {min(x['s'] for x in r)*G:.2f} to {max(x['s'] for x in r)*G:.2f} Gpc\nall seen light was emitted before our crossing ({m['tc']*Gy:.2f} Gyr)",fontsize=11,loc='left')
h=[plt.Line2D([],[],color='blue',lw=5),plt.Line2D([],[],color='r',lw=3),plt.Line2D([],[],color='g',lw=2),plt.Line2D([],[],color='darkorange'),plt.Line2D([],[],color='darkorange',ls=':'),plt.Line2D([],[],color='gray',ls='--'),plt.Line2D([],[],color='m',marker='.',ls='')]
fig.legend(h,['blue arc: emission points on the recombination sphere (the only part we can see)','red: same emitters today','green: c × time since we crossed, centred on us today','orange: straight light path, emission point → us today (every 2.5° of sky)','orange dotted: path whose redshift comes out negative or infinite','gray: dashed radial line from CoE, arrow = emitter moving out to today','magenta: emission points where 1+z ≤ 0 or infinite'],loc='lower left',bbox_to_anchor=(0.02,1-1.6/FH),fontsize=11,frameon=False,ncol=2)
fig.suptitle('Static space (no expansion): matter moves through fixed space, light travels straight at c; recombination radius 0.5 D; D set so green = 4.2 Gpc',y=1-0.25/FH,fontsize=13)
fig.savefig('static_rays_2x2.png',dpi=90)
print('axes sizes',set(tuple(np.round(a.get_position().bounds[2:],4)) for a in axs.ravel()))
