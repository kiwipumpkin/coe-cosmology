# Mean direction of the published dipole ends (6 observed pairs, CoE model excluded), smallest covering circle, and the antipodal circle; drawn on the Mollweide.
import os,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
H=os.path.dirname(os.path.abspath(__file__))
src=open(H+'/dipoles_sky_map_lines.py').read();exec(src.split('def xy')[0].split('\n',3)[3])  # P
def v(l,b):l,b=np.radians(l),np.radians(b);return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def lb(u):u=u/np.linalg.norm(u);return np.degrees(np.arctan2(u[1],u[0]))%360,np.degrees(np.arcsin(u[2]))
ang=lambda a,b:np.degrees(np.arccos(np.clip(a@b,-1,1)))
obs=[p for p in P if p[5] and 'CoE' not in p[0]]
V=np.array([v(p[1],p[2]) for p in obs])
m=V.mean(0);m/=np.linalg.norm(m)
# minimal enclosing circle on sphere: refine centre by minimising max angle
best=m;rb=max(ang(m,x) for x in V)
for it in range(20000):
  c=best+np.random.default_rng(it).normal(0,.03,3);c/=np.linalg.norm(c);r=max(ang(c,x) for x in V)
  if r<rb:best,rb=c,r
lm,bm=lb(m);lc,bc=lb(best)
rmean=max(ang(m,x) for x in V)
coe=v(180,-15)
res=dict(members=[p[0] for p in obs],mean_dir=[round(lm,1),round(bm,1)],mean_dir_cover_radius=round(rmean,1),
 min_circle_centre=[round(lc,1),round(bc,1)],min_circle_radius=round(rb,1),antipode_centre=[round((lc+180)%360,1),round(-bc,1)],
 seps_from_centre={p[0]:round(ang(best,x),1) for p,x in zip(obs,V)},CoE_sep_from_centre=round(ang(best,coe),1),AoE_sep=round(ang(best,v(260,60)),1))
json.dump(res,open(H+'/dipole_cluster_circle.json','w'),indent=1);print(json.dumps(res,indent=1))
def xy(l,b):return -np.radians(((np.asarray(l)+180)%360)-180),np.radians(b)
def circ(c,r):
  c=c/np.linalg.norm(c);a=np.cross(c,[0,0,1.]);a/=np.linalg.norm(a);b2=np.cross(c,a);t=np.linspace(0,2*np.pi,721)
  pts=np.cos(np.radians(r))*c[:,None]+np.sin(np.radians(r))*(np.outer(a,np.cos(t))+np.outer(b2,np.sin(t)))
  L=np.degrees(np.arctan2(pts[1],pts[0]))%360;B=np.degrees(np.arcsin(pts[2]));return xy(L,B)
fig=plt.figure(figsize=(12,6.6));ax=fig.add_subplot(111,projection='mollweide');ax.grid(alpha=.4)
for n,l,b,c,mk,anti in P:
  x,y=xy(l,b);ax.plot(x,y,mk,color=c,mec='k',ms=14 if mk!='*' else 20,label=f'{n} ({l:g}°, {b:+g}°)')
  if anti:x2,y2=xy(l+180,-b);ax.plot([x,x2],[y,y2],'-',color=c,lw=2,alpha=.8);ax.plot(x2,y2,mk,mfc='none',mec=c,mew=2,ms=12 if mk!='*' else 18)
for cc,col,lab in [(best,'k','published ends'),(-best,'gray','antipodes')]:
  X,Y=circ(cc,rb);br=np.where(np.abs(np.diff(X))>np.pi)[0]+1
  for i,(xs,ys) in enumerate(zip(np.split(X,br),np.split(Y,br))):ax.plot(xs,ys,'--',color=col,lw=2.5,label=f'{rb:.1f}° circle, {lab}' if i==0 else None)
  ax.plot(*xy(*lb(cc)),'+',color=col,ms=22,mew=3)
ax.set_xticks(np.radians([-150,-120,-90,-60,-30,0,30,60,90,120,150]));ax.set_xticklabels(['150','120','90','60','30','0','330','300','270','240','210'])
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.05),ncol=3,fontsize=9,frameon=False)
ax.set_title(f'Smallest circle covering the 6 observed published ends: centre (l,b)=({lc:.1f}°, {bc:+.1f}°), radius {rb:.1f}°; antipodal circle dashed grey',fontsize=10,pad=14)
fig.tight_layout();fig.savefig(H+'/dipole_cluster_circle.png',dpi=140,facecolor='white',bbox_inches='tight')
