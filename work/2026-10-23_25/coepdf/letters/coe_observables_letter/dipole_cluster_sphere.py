# Sphere views of the dipole pairs with the 44.5° covering circle (published ends) and its antipodal circle, from dipole_cluster_circle.json.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
src=open(H+'/dipoles_sky_map_lines.py').read();exec(src.split('def xy')[0].split('\n',3)[3])
R=json.load(open(H+'/dipole_cluster_circle.json'))
def v(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
C=v(*R['min_circle_centre']);r=np.radians(R['min_circle_radius'])
a=np.cross(C,[0,0,1.]);a/=np.linalg.norm(a);b2=np.cross(C,a);t=np.linspace(0,2*np.pi,361)
ring=lambda c:np.cos(r)*c[:,None]+np.sin(r)*(np.outer(a,np.cos(t))+np.outer(b2,np.sin(t)))
views=[(0,R['min_circle_centre'][0]+90,'side-on to the axis'),(R['min_circle_centre'][1],R['min_circle_centre'][0],'looking straight down the axis (+ end)'),(20,200,'az 200, el 20'),(20,20,'az 20, el 20'),(60,240,'from above, el 60'),(-60,60,'from below, el -60')]
fig=plt.figure(figsize=(18,12))
for k,(el,az,ttl) in enumerate(views):
  ax=fig.add_subplot(2,3,k+1,projection='3d')
  for B in (-60,-30,0,30,60):
    rr=np.cos(np.radians(B));ax.plot(rr*np.cos(t),rr*np.sin(t),np.sin(np.radians(B))+0*t,'k',lw=1.1 if B==0 else .4,alpha=.6 if B==0 else .25)
  for L in range(0,360,30):
    bb=np.radians(np.linspace(-90,90,60));ax.plot(np.cos(bb)*np.cos(np.radians(L)),np.cos(bb)*np.sin(np.radians(L)),np.sin(bb),'k',lw=.4,alpha=.2)
  for L,lab in [(0,'l=0'),(90,'l=90'),(180,'l=180'),(270,'l=270')]:ax.text(*(1.15*v(L,0)),lab,fontsize=8,color='gray')
  ax.text(0,0,1.12,'NGP',fontsize=8,color='gray')
  ax.plot(*ring(C),'k--',lw=2.5);ax.plot(*ring(-C),'--',color='gray',lw=2.5)
  ax.plot([C[0],-C[0]],[C[1],-C[1]],[C[2],-C[2]],color='k',lw=1,alpha=.5);ax.scatter(*C,marker='+',s=300,color='k',linewidths=3);ax.scatter(*-C,marker='+',s=300,color='gray',linewidths=3)
  for n,l,b,c,mk,anti in P:
    p=v(l,b);ax.scatter(*p,color=c,marker=mk,s=120 if mk!='*' else 260,edgecolors='k',depthshade=False,label=n if k==0 else None)
    if anti:ax.plot(*np.array([p,-p]).T,color=c,lw=1.8);ax.scatter(*-p,facecolors='none',edgecolors=c,marker=mk,s=110 if mk!='*' else 230,linewidths=2,depthshade=False)
  ax.set_box_aspect((1,1,1));[f(-1,1) for f in (ax.set_xlim,ax.set_ylim,ax.set_zlim)];ax.view_init(el,az);ax.set_axis_off();ax.set_title(ttl,fontsize=11)
fig.legend(loc='lower center',ncol=5,fontsize=10,frameon=False)
fig.suptitle(f"Covering circle (black dashed) centre ({R['min_circle_centre'][0]}°, {R['min_circle_centre'][1]:+}°), radius {R['min_circle_radius']}°; antipodal circle (grey dashed) at ({R['antipode_centre'][0]}°, {R['antipode_centre'][1]:+}°)",fontsize=13)
fig.savefig(H+'/dipole_cluster_sphere.png',dpi=110,facecolor='white',bbox_inches='tight')
