import pickle,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':15})
P=pickle.load(open('rays_exp.pkl','rb'));C=json.load(open('res_exp.json'))['cmb_A']
fig,ax=plt.subplots(1,2,figsize=(18,9))
a=ax[0];a.set_aspect('equal')
a.scatter(0,0,marker='x',s=300,c='r');a.text(0.1,-0.6,'CoE',color='r')
for R in (3,6):a.add_patch(plt.Circle((0,0),R,fill=False,ls=':',color='0.5'))
a.add_patch(plt.Circle((0,0),11,fill=False,color='g',lw=2));a.text(-10.5,8.5,'outer edge of matter',color='g',fontsize=13)
ox=8;a.scatter(ox,0,marker='*',s=500,c='k');a.text(ox+0.3,0.4,'us (now)')
for t in np.linspace(0,2*np.pi,13)[:-1]:
    d=np.array([np.cos(t),np.sin(t)]);o=np.array([ox,0.])
    bq=o@d;L=-bq+np.sqrt(bq*bq-(o@o-121))
    e=o+L*d;a.plot([ox,e[0]],[0,e[1]],color='0.6')
    bad=e[0]<0 and abs(e[1])<11
    a.scatter(*e,s=90,c='r' if bad else 'b',zorder=3)
a.scatter([],[],c='b',label='emitter at matter edge');a.scatter([],[],c='r',label='emitter past the CoE (forbidden)')
a.legend(loc='lower left',fontsize=12);a.set_xlim(-12,12);a.set_ylim(-12,12)
a.set_title('What the code did\nNo element reaches 3000 K, so each ray runs out\nto the edge of the matter and that point emits',fontsize=14)
b=ax[1];b.set_aspect('equal')
b.scatter(0,0,marker='x',s=300,c='r');b.text(0.1,-0.5,'CoE',color='r')
for R in (3,6):b.add_patch(plt.Circle((0,0),R,fill=False,ls=':',color='0.5'))
ox=8;b.scatter(ox,0,marker='*',s=500,c='k');b.text(ox+0.3,0.4,'us (now)')
b.add_patch(plt.Circle((ox,0),2.5,fill=False,color='orange',lw=3))
for t in np.linspace(0,2*np.pi,13)[:-1]:
    if np.cos(t)<-0.6: continue
    b.annotate('',xy=(ox,0),xytext=(ox+2.5*np.cos(t),2.5*np.sin(t)),arrowprops=dict(arrowstyle='->',color='orange'))
b.text(ox-3.8,2.9,'emission circle C(n-1):\nlight leaves at t>=t*, all on\nour past light cone',fontsize=13,color='darkorange')
b.text(-5,-6,'CoE side cut off: horizon',fontsize=13,color='r')
b.set_xlim(-7,12);b.set_ylim(-7,7)
b.set_title('What we drew\nEmitters sit on the circle around us,\nnone reach past the CoE',fontsize=14)
fig.text(0.5,0.02,f"Cause: K = {C['K']:.2f} K, so no element ever reaches 3000 K on the ray and the t* rule never fires. The code falls back to the matter edge. 11/60 emitters lie past the CoE.",ha='center',fontsize=14)
plt.savefig('bug_sketch.png',dpi=90,bbox_inches='tight')
