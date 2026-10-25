import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;dr=0.07;N=6
fig,ax=plt.subplots(figsize=(16,10));ax.set_aspect('equal');ax.axis('off')
ax.add_patch(plt.Circle((0,0),Rrec,color='#ffd6d6',zorder=0))
ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
ax.plot(0,0,'kx',ms=14,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=16)
cols=plt.cm.viridis(np.linspace(0,0.9,N))
for i in range(N):
    xu=Rrec+U*dr*i
    ax.plot(xu,0,'*',color=cols[i],ms=20,zorder=5)
    ax.text(xu-0.02,-0.05,str(i),fontsize=14,color=cols[i])
    if i>0: ax.add_patch(plt.Circle((xu,0),dr*i,fill=False,color=cols[i],lw=2.5))
ax.set_xlim(-0.2,2.0);ax.set_ylim(-0.5,0.5)
ax.set_title('SCHEMATIC: all steps. Star n = us at step n (U=1.5c away from CoE); circle n = visible radius c*t_n around us\nStep 0 is on the recombination line (radius 0). At step n we see light emitted from circle n-1',fontsize=16)
plt.savefig('all_steps.png',dpi=100,bbox_inches='tight')
