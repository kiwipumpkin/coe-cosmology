import sys,io,contextlib;sys.path.insert(0,'../..');sys.path.insert(0,'..')
with contextlib.redirect_stdout(io.StringIO()):
    import fit,free
import numpy as np,json
from scipy.spatial.transform import Rotation as R
vec=fit.vec;lb=fit.lb
S0=json.load(open('../../rotfit_stats.json'))
e1,e2=vec(0,15),vec(90,0);c=np.radians(np.arange(0,360,.25));N=np.outer(np.cos(c),e1)+np.outer(np.sin(c),e2)
L=np.degrees(np.arctan2(N[:,1],N[:,0]))%360;o=np.argsort(L);L=L[o];N=N[o]
a0=vec(0,15);v0=vec(*S0['fit']['motion_lb'])
def curve(a,v):
    Vs=fit.V;fit.V=N;m=fit.sky(a,v);fit.V=Vs;return m-free.c0-free.c1*(N@a)
def dl(x,y):return (x-y+180)%360-180
def pk(y):return L[np.argmax(y)],L[np.argmin(y)]
def rot(seq,angs):return R.from_euler(seq,angs,degrees=True).as_matrix()
h0,c0=pk(curve(a0,v0))
res={'base':dict(hot=float(h0),cold=float(c0))}
sw={}
for ax in 'xyz':
    for t in list(range(-90,91,5))+list(range(-10,11,1)):
        Rm=rot(ax,t);h,cl=pk(curve(Rm@a0,Rm@v0));sw[f'{ax}{t}']=dict(hot=float(h),cold=float(cl),dhot=float(dl(h,h0)),dcold=float(dl(cl,c0)))
res['single']=sw
# 2D grid search: rotations about any axis through observer = all orientations; search yaw(z) then tilt x then tilt y
best=[];
for z in range(0,360,15):
  for x in range(-90,91,15):
    for y in range(-90,91,15):
      Rm=rot('zxy',[z,x,y]);h,cl=pk(curve(Rm@a0,Rm@v0));best.append((max(abs(dl(h,264.02)),abs(dl(cl,84.02))),z,x,y,h,cl))
best.sort();res['grid_top']=[list(map(float,b)) for b in best[:10]]
print(json.dumps({k:v for k,v in res.items() if k!='single'},indent=0))
for k,v in sw.items():
    if k[0]!='z' or True: pass
json.dump(res,open('spread.json','w'),indent=1)
