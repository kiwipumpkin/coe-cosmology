import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;c=1.0
fig,axs=plt.subplots(2,3,figsize=(18,12))
th=np.linspace(0,2*np.pi,400)
for i,ax in enumerate(axs.flat):
    dt=0.07*i
    ax.set_aspect('equal');ax.axis('off')
    ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
    ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
    ax.plot(0,0,'kx',ms=12,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=13)
    x0=Rrec;xu=Rrec+U*dt
    ax.plot([x0,xu],[0,0],'b:',lw=2)
    ax.plot(x0,0,'o',color='gray',ms=8);ax.text(x0-.05,0.05,'t*',fontsize=13,color='gray')
    ax.plot(xu,0,'b*',ms=20)
    if dt>0: ax.annotate('',xy=(xu+0.15,0),xytext=(xu,0),arrowprops=dict(arrowstyle='->',color='blue',lw=3))
    r=c*dt
    ax.plot(x0+r*np.cos(th),r*np.sin(th),color='orange',lw=3)
    ax.set_title(f'Panel {i+1}: t-t* = {dt:.2f}\nus at U={U}c, light front at c',fontsize=15)
    ax.set_xlim(-0.2,1.9);ax.set_ylim(-0.6,0.6)
fig.suptitle('SCHEMATIC: light emitted at t* where we crossed recombination (orange, grows at c)\nwhile we (blue star) move away from the CoE at U=1.5c, so we outrun the front on the CoE side',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.9]);plt.savefig('six_motion.png',dpi=100)
