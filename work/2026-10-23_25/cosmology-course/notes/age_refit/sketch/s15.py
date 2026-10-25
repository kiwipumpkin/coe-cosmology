import matplotlib.pyplot as plt, numpy as np
Rrec=0.65;U=1.5;dr=0.07;N=7
cols=['red','orange','gold','green','blue','indigo','violet']
fig,ax=plt.subplots(figsize=(18,10));ax.set_aspect('equal');ax.axis('off')
ax.add_patch(plt.Circle((0,0),Rrec,color='#ffe8e8',zorder=0))
ax.add_patch(plt.Circle((0,0),Rrec,fill=False,color='crimson',lw=2,ls='--'))
ax.plot(0,0,'kx',ms=14,mew=2);ax.text(0.02,-0.08,'CoE',fontsize=16)
for i in range(N):
    xu=Rrec+U*dr*i;r=dr*i
    if i>0:
        ax.add_patch(plt.Circle((xu,0),r,fill=False,color=cols[i],lw=2.5))
        ax.text(xu,r+0.012,f'circle {i}',color=cols[i],fontsize=13,ha='center',fontweight='bold')
    oc='black' if i==0 else cols[i-1]
    ax.plot(xu,0,'*',color=oc,ms=22,mec='k',mew=0.8,zorder=6)
    ax.text(xu,-0.06-0.012*(i%2),f'obs {i}',color=oc,fontsize=13,ha='center',fontweight='bold')
ax.set_xlim(-0.2,2.3);ax.set_ylim(-0.55,0.55)
ax.set_title('SCHEMATIC: circle n (radius c*t_n around observer n). Observer n is coloured like circle n-1,\nthe light it sees. Obs 0 (black) sits on the recombination line and sees nothing. U=1.5c away from CoE',fontsize=16)
plt.savefig('all_steps_numbered.png',dpi=100,bbox_inches='tight')
