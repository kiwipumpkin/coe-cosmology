import numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':18})
cols=['red','orange','gold','green','blue','indigo','violet']
fig,a=plt.subplots(figsize=(11,11));a.set_aspect('equal')
for i,c in enumerate(cols):
    t=i+1;Ri,Ro=1.5*t,2.0*t
    th=np.linspace(0,2*np.pi,400)
    a.fill(np.r_[Ro*np.cos(th),Ri*np.cos(th[::-1])],np.r_[Ro*np.sin(th),Ri*np.sin(th[::-1])],color=c,alpha=0.12)
    a.plot(Ro*np.cos(th),Ro*np.sin(th),color=c,lw=2);a.plot(Ri*np.cos(th),Ri*np.sin(th),color=c,lw=1,ls='--')
    a.text(Ro*0.71+0.1,Ro*0.71,f't{t}',color=c,fontsize=16)
for i,c in enumerate(cols):
    t=i+1;r=1.75*t;a.scatter(r,0,s=160,c=c,edgecolors="k",zorder=5)
a.text(1.75*7+0.3,0.5,"us",fontsize=16)
a.scatter(0,0,marker='x',s=400,c='k');a.text(0.3,-0.6,'CoE')
for ang in np.radians(range(0,360,45)):a.annotate('',xy=(13*np.cos(ang),13*np.sin(ang)),xytext=(10*np.cos(ang),10*np.sin(ang)),arrowprops=dict(arrowstyle='->',lw=2))
a.set_xlim(-15,15);a.set_ylim(-15,15);a.set_xticks([]);a.set_yticks([])
a.set_title('Our shell expanding outward from the CoE\n(one colour per time step; arrows: shell moves out in every direction)')
plt.savefig('shell1.png',dpi=80,bbox_inches='tight')
