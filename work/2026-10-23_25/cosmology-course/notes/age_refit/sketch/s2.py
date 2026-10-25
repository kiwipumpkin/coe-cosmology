import matplotlib.pyplot as plt, numpy as np
fig,ax=plt.subplots(figsize=(16,10));ax.set_aspect('equal')
D=1.0;Rrec=D-0.35
ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=3,ls='--'))
for r in [D,D+0.35]:
    ax.add_patch(plt.Circle((0,0),r,fill=False,color='gray',lw=3 if r==D else 2,ls='-' if r==D else '--'))
ax.text(-0.3,0.25,'PLASMA',fontsize=22,color='crimson')
ax.text(-0.5,-Rrec-0.12,'recombination surface',fontsize=18,color='crimson')
ax.plot(0,0,'kx',ms=18,mew=3);ax.text(0.03,-0.1,'CoE',fontsize=20)
ax.plot(D,0,'b*',ms=26);ax.text(D+.03,-.1,'us',fontsize=20,color='b')
th=np.linspace(0,2*np.pi,400)
for r in [0.1,0.2,0.3]:
    x,y=D+r*np.cos(th),r*np.sin(th);m=np.hypot(x,y)>=Rrec
    ax.plot(np.where(m,x,np.nan),np.where(m,y,np.nan),color='orange',lw=3)
    ax.plot(np.where(~m,x,np.nan),np.where(~m,y,np.nan),color='orange',lw=1,ls=':')
for t in np.linspace(0,2*np.pi,12,endpoint=False):
    x,y=D+0.3*np.cos(t),0.3*np.sin(t);u=np.array([x,y])/np.hypot(x,y)
    ax.arrow(x,y,.1*u[0],.1*u[1],head_width=.025,color='red')
ax.text(-1.55,1.45,'Gray = SHELL (around CoE, expanding outward)\nRed dashed = recombination surface (inside = plasma)\nOrange = EMISSION CIRCLES around us, c(t-t*); solid = recombined, dotted = still plasma\nRed arrows = emitter motion away from CoE',fontsize=17,va='top')
ax.set_xlim(-1.6,1.8);ax.set_ylim(-1.5,1.5);ax.axis('off');ax.set_title('SCHEMATIC (not computed)',fontsize=22)
plt.savefig('geom2.png',dpi=100,bbox_inches='tight')
