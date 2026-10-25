import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;dr=0.07;N=6
cols=plt.cm.plasma(np.linspace(0,0.85,N))
fig,axs=plt.subplots(N,1,figsize=(12,36))
angs=np.radians([0,45,90,135,180,225,270,315])
for i,ax in enumerate(axs):
    ax.set_aspect('equal');ax.axis('off')
    ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe0e0',zorder=0))
    ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
    ax.plot(0,0,'kx',ms=12,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=13)
    xu=Rrec+U*dr*i
    ax.plot(xu,0,'*',color='blue',ms=20,zorder=6)
    if i>0: ax.add_patch(plt.Circle((xu,0),dr*i,fill=False,color=cols[i],lw=3))
    if i>1:
        xp=Rrec+U*dr*(i-1);r=dr*(i-1)
        ax.add_patch(plt.Circle((xp,0),r,fill=False,color=cols[i-1],lw=2.5,ls='--'))
        for a in angs:
            ex,ey=xp+r*np.cos(a),r*np.sin(a);ax.plot(ex,ey,'o',color=cols[i-1],ms=7)
            ax.annotate('',xy=(xu,0),xytext=(ex,ey),arrowprops=dict(arrowstyle='->',color=cols[i-1],lw=2))
    t='Step 0: on the recombination line, visible distance c*0, no light' if i==0 else (f'Step 1: circle 1 (solid); light from step 0 = a point, nothing yet' if i==1 else f'Step {i}: circle {i} (solid); light arriving in step-{i-1} colour from circle {i-1} (dashed)')
    ax.set_title(t,fontsize=15,color=cols[i] if i>0 else 'k')
    ax.set_xlim(-0.2,1.9);ax.set_ylim(-0.5,0.5)
fig.suptitle('SCHEMATIC: each step has its own colour (time). Light reaching us carries the previous step\'s colour. U=1.5c away from CoE',fontsize=16)
plt.tight_layout(rect=[0,0,1,0.98]);plt.savefig('six_v9.png',dpi=90)
