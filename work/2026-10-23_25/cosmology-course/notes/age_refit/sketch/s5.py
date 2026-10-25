import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;dr=0.07
fig,axs=plt.subplots(2,3,figsize=(18,12))
angs=np.radians([0,45,90,135,180,225,270,315])
for i,ax in enumerate(axs.flat):
    ax.set_aspect('equal');ax.axis('off')
    ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
    ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
    ax.plot(0,0,'kx',ms=12,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=13)
    xu=Rrec+U*dr*i
    ax.plot(xu,0,'b*',ms=20,zorder=5)
    ax.annotate('',xy=(xu+0.2,0.0),xytext=(xu+0.02,0),arrowprops=dict(arrowstyle='->',color='blue',lw=3))
    for k in range(1,i+1):
        r=dr*k
        ax.add_patch(plt.Circle((xu,0),r,fill=False,color='orange',lw=3 if k==i else 1.2,alpha=1 if k==i else .6))
    if i>0:
        r=dr*i
        for a in angs:
            ax.annotate('',xy=(xu+0.02*np.cos(a),0.02*np.sin(a)),xytext=(xu+r*np.cos(a),r*np.sin(a)),arrowprops=dict(arrowstyle='->',color='green',lw=1.8))
    ax.set_title(f'Step {i}: we are {xu:.2f} from CoE\nseeing circle {i} (radius {dr*i:.2f})',fontsize=15)
    ax.set_xlim(-0.2,1.9);ax.set_ylim(-0.6,0.6)
fig.suptitle('SCHEMATIC: emission circles centred on us, travelling with us (U=1.5c away from CoE)\nEach step we are farther out and the next circle (bold) becomes visible; green = light reaching us',fontsize=17)
plt.tight_layout(rect=[0,0,1,0.9]);plt.savefig('six_v3.png',dpi=100)
