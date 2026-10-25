import matplotlib.pyplot as plt, numpy as np
fig,ax=plt.subplots(figsize=(16,10));ax.set_aspect('equal')
D=1.0
for r,a in [(D-0.35,.5),(D,1),(D+0.35,.5)]:
    ax.add_patch(plt.Circle((0,0),r,fill=False,color='gray',lw=2*a+1,ls='-' if a==1 else '--'))
ax.plot(0,0,'kx',ms=18,mew=3);ax.text(0.03,-0.1,'CoE',fontsize=20)
ox,oy=D,0;ax.plot(ox,oy,'b*',ms=26);ax.text(ox+.03,oy-.1,'us',fontsize=20,color='b')
for r in [0.1,0.2,0.3]:
    ax.add_patch(plt.Circle((ox,oy),r,fill=False,color='orange',lw=2.5))
for th in np.linspace(0,2*np.pi,12,endpoint=False):
    x,y=ox+0.3*np.cos(th),oy+0.3*np.sin(th);u=np.array([x,y])/np.hypot(x,y)
    ax.arrow(x,y,.1*u[0],.1*u[1],head_width=.025,color='red')
ax.text(-1.55,1.45,'Gray = SHELL (concentric around CoE, expanding outward)\nOrange = EMISSION CIRCLES (concentric around us, radius c(t-t*))\nRed arrows = emitter motion away from CoE',fontsize=18,va='top')
ax.set_xlim(-1.6,1.8);ax.set_ylim(-1.5,1.5);ax.axis('off')
ax.set_title('SCHEMATIC (not computed)',fontsize=22)
plt.savefig('geom.png',dpi=100,bbox_inches='tight')
