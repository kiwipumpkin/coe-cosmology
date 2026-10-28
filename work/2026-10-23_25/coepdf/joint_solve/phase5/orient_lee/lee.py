# Look-elsewhere MC for orientation of CoE model. Reuses 23_sky_mapA/B_al.py conventions exactly:
# residual = interp(model-frame longitude, M25['l'], M25[A|B], period=360); match = (val>0)==(hot).
# Paper orientation R = rot(coe,45) @ rot(y,-15); model frame m = v @ R ; sky = R @ m.
import json,sys,numpy as np
H=__file__.rsplit('/',1)[0]
MM=json.load(open(H+'/M25.json'))
exec(open('/workspace/coepdf/book/code/24_spots.py').read())
SP=[(K,)+s for K in 'AB' for s in SPOTS[K]]
def vec(l,b):l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
def rot(ax,a):
  ax=ax/np.linalg.norm(ax);a=np.radians(a);K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]])
  return np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K
R1=rot(np.array([0,1.,0]),-15.0);coe=R1@vec(180,0);R0=rot(coe,45.0)@R1
V=np.array([vec(s[2],s[3]) for s in SP]);HOT=np.array([s[4]=='hot' for s in SP]);KA=np.array([s[0]=='A' for s in SP])
NAMES=[s[1] for s in SP]
PAIRS=[('CMB dipole cold','CMB dipole hot'),('Bulk flow away','Bulk flow toward'),('NVSS cold','NVSS hot'),('Hemi. asym. cold','Hemi. asym. hot'),('CatWISE cold','CatWISE hot')]
PI=[(NAMES.index(a),NAMES.index(b)) for a,b in PAIRS];SINGLE=[NAMES.index('CMB Cold Spot'),NAMES.index('Loop I / NPS')]
AOE=vec(260.,60.)
def matches(Rs):           # Rs (N,3,3) -> (N,12) bool
  m=np.einsum('sj,njk->nsk',V,Rs);lm=np.degrees(np.arctan2(m[...,1],m[...,0]))%360
  val=np.where(KA,np.interp(lm,MM['l'],MM['A'],period=360),np.interp(lm,MM['l'],MM['B'],period=360))
  return (val>0)==HOT,val
def paired(M): return np.stack([M[:,a]&M[:,b] for a,b in PI]+[M[:,i] for i in SINGLE],1).sum(1)
def aoe_off(Rs,DR):        # critical pts: model-frame b=0, l=180±theta, theta=arccos(R/(2D)); plus antipodes
  th=np.arccos(1/(2*DR));P=vec(np.array([180+np.degrees(th),180-np.degrees(th)]),np.zeros(2))
  s=np.einsum('njk,pk->npj',Rs,P)@AOE;return np.degrees(np.arccos(np.clip(np.abs(s).max(1),-1,1)))
def haar(n,rng):
  q=rng.normal(size=(n,4));q/=np.linalg.norm(q,axis=1)[:,None];w,x,y,z=q.T
  return np.stack([1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w),2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w),2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)],1).reshape(-1,3,3)
if __name__=='__main__':
  N=int(sys.argv[1]) if len(sys.argv)>1 else 200000;rng=np.random.default_rng(12345)
  DRS={'phase4b joint D=27.34/R_ls=21.74':27.34/21.74,'phase4 joint D=27.92/21.77':27.92/21.77,'D=42/21.74':42/21.74,'D=58.2/21.74':58.2/21.74,'figure D/R=5':5.0}
  M0,v0=matches(R0[None]);s0=int(M0.sum());p0=int(paired(M0)[0])
  out={'real':{'score12':s0,'scoreA':int(M0[0,KA].sum()),'scoreB':int(M0[0,~KA].sum()),'paired7':p0,'vals':dict(zip(NAMES,np.round(v0[0],3).tolist()))},'N':N}
  Rh=haar(N,rng);phi=rng.uniform(0,360,N);Rr=np.array([rot(coe,p)@R1 for p in phi])
  res={}
  for tag,Rs in [('isotropic',Rh),('fixedCoE_roll',Rr)]:
    M,_=matches(Rs);S=M.sum(1);P=paired(M);r={'score_hist':np.bincount(S,minlength=13).tolist(),'p_score>=real':float((S>=s0).mean()),'p_paired>=real':float((P>=p0).mean()),'paired_hist':np.bincount(P,minlength=8).tolist(),'aoe':{}}
    for dn,DR in DRS.items():
      o0=float(aoe_off(R0[None],DR)[0]);o=aoe_off(Rs,DR)
      r['aoe'][dn]={'DR':DR,'real_offset_deg':o0,'p_aoe<=real':float((o<=o0).mean()),'p_joint_score':float(((S>=s0)&(o<=o0)).mean()),'p_joint_paired':float(((P>=p0)&(o<=o0)).mean())}
      if dn.startswith('phase4b'):np.save(H+f'/off_{tag}.npy',o.astype(np.float32))
    np.save(H+f'/score_{tag}.npy',S.astype(np.int8));res[tag]=r
  out['mc']=res;json.dump(out,open(H+'/lee_results.json','w'),indent=1);print(json.dumps(out,indent=1))
