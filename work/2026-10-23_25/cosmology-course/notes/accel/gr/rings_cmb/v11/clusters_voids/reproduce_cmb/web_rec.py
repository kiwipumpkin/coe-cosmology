import os,sys,json,numpy as np,healpy as hp
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree
sys.path.insert(0,'../../rotfit');import fit
SEQ='/workspace/clump_scan/sequential';c=os.getcwd();os.chdir(SEQ);sys.path.insert(0,SEQ);import seq_core as S;os.chdir(c)
DGPC=float(sys.argv[1]) if len(sys.argv)>1 else 42.0; TAG=sys.argv[2] if len(sys.argv)>2 else 'D42'
nside=256;V=np.array(hp.pix2vec(nside,np.arange(hp.nside2npix(nside)))).T
o=np.load('../../shrink_rec/b_1100_k0.npz')['o']   # psi, 1+z, t_e, lookback, w (emitter label / D), phi (angle at CoE), R
aG=fit.vec(0,15)                                    # anti-CoE, galactic
psi=np.degrees(np.arccos(np.clip(V@aG,-1,1)))
smooth=np.exp(-CubicSpline(o[:,0],np.log(o[:,1]))(psi));smooth=smooth/smooth.mean()-1
wf=CubicSpline(o[:,0],o[:,4]);pf=CubicSpline(o[:,0],o[:,5])
# web: Chris's frame -> galactic. Rotation taking his anti-CoE (-C_G) to our anti-CoE (0,+15); roll about the axis is arbitrary (minimal rotation)
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X,r,M=f['Xp'],f['rp'],f['Mp'];isV=f['kp']=='void'
CG=np.array(S.C_G);ah=-CG/np.linalg.norm(CG)
from scipy.spatial.transform import Rotation as Rr
k=np.cross(ah,aG);s=np.linalg.norm(k);ang=np.arctan2(s,ah@aG);Rm=Rr.from_rotvec(k/s*ang).as_matrix()
Xg=X@Rm.T
rho=S.rho_gpc(X);GPC=3.0857e25
dlt=M/(4/3*np.pi*(r*GPC)**3)/rho if np.nanmax(rho)<1e-10 else M/(4/3*np.pi*r**3)/rho   # signed contrast (voids negative)
# --- grow back to recombination: nonlinear -> linear (Bernardeau 1994 spherical approx), then linear growth D∝a
ZM=float(np.mean(o[:,1]));dNL=dlt.copy();dlt=1.5*(1-(1+dNL)**(-2/3))/ZM
print('delta_rec void',np.percentile(dlt[isV],[0,50,100]),'cluster',np.percentile(dlt[~isV],[0,50,100]))
# emission point (today's comoving position) for each sky direction: E = -D a + w D (cos phi a + sin phi t)
t=V-(V@aG)[:,None]*aG;tn=np.linalg.norm(t,axis=1);t=np.where(tn[:,None]>1e-9,t/np.maximum(tn,1e-12)[:,None],0)
w=wf(psi);ph=pf(psi);E=DGPC*(-aG[None,:]+w[:,None]*(np.cos(ph)[:,None]*aG[None,:]+np.sin(ph)[:,None]*t))
dist=np.linalg.norm(E,axis=1)
tree=cKDTree(Xg);hits=tree.query_ball_point(E,r.max())
web=np.zeros(len(V));nhit=np.zeros(len(V),int);used=set()
for i,h in enumerate(hits):
    if not h:continue
    h=np.array(h);d=np.linalg.norm(Xg[h]-E[i],axis=1);h=h[d<r[h]]
    if len(h):web[i]=np.sum(np.log1p(dlt[h]))/3;nhit[i]=len(h);used.update(h.tolist())
print("delta void",np.percentile(dlt[isV],[0,50,100]),"cluster",np.percentile(dlt[~isV],[0,50,100]))
tot=smooth+web
def strip(x):return hp.remove_dipole(x)
def lrms(x):cl=hp.anafast(x,lmax=12);l=np.arange(13);return np.sqrt((2*l+1)*cl/(4*np.pi))/1e-5
out={'D_Gpc':DGPC,'surface_distance_from_us_Gpc':[float(dist.min()),float(np.median(dist)),float(dist.max())],
 'bodies_touching_surface':len(used),'voids_touching':int(sum(isV[list(used)])),'clusters_touching':int(sum(~isV[list(used)])),
 'sky_fraction_covered':float(np.mean(nhit>0)),'sky_fraction_in_void':float(np.mean(web<0)),'sky_fraction_in_cluster':float(np.mean(web>0))}
for name,x in (('smooth',smooth),('web_only',web),('smooth_plus_web',tot)):
    y=strip(x);L=lrms(y);g=lambda ls:float(np.sqrt(np.sum(L[ls]**2)))
    out[name]=dict(map_rms_1e5=float(np.std(y)/1e-5),map_rms_uK=float(np.std(y)*2.7255e6),dipole_before_motion_xPlanck=float(np.linalg.norm(hp.fit_dipole(x)[1])/(3362/2.7255e6)),
      l2=float(L[2]),l3_5=g([3,4,5]),l6_12=g(list(range(6,13))),per_l=[float(v) for v in L[1:]])
print(json.dumps(out,indent=1));json.dump(out,open(f'web_rec_{TAG}.json','w'),indent=1)
np.save(f"finalmap_{TAG}.npy",np.asarray(strip(tot)));np.save(f'webrec_{TAG}.npy',web)
