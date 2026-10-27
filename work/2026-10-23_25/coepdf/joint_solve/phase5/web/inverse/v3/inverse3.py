# inverse3.py (v3; copied/derived from ../inverse2.py, which is unchanged). Chris, Oct 26: start from the SAME statistically uniform random
# 3000-body web of clusters+voids as v2 (identical rng draw order -> identical y0, delta0, candidates for a given seed), on the Route A
# eps=0.01 St background, and solve the inverse problem with:
#  (1) objective = weighted alm residual (bg+web) - Planck target, l=1..LC with LC=130 (v2: 40 fitted + 41..60 penalised to 0).
#      Target = SMICA l=2..130 at Nside-64 pixel centres (make_target.py; l=1 target 0 as in v2).
#  (2) each body's density contrast delta (mass M = delta*rho_c*V) free alongside its position. Only physical bound: delta >= -1
#      (empty void). NO upper cap on delta. Geometric limits only: r in [RMIN,RMAXF]*R, RMAXF=0.975 = one softening length inside the
#      CMB shell (kernel softening EPS=0.025 R), RMIN=0.01 avoids r=0. Bodies at any of these limits are counted and flagged.
#  (3) sparsity on the NUMBER of bodies changed (L0, greedy): each greedy step changes ONE body, choosing among (a) move to one of NC random
#      candidate positions with its optimal delta (closed form, clipped at -1) and (b) keep position, re-optimise delta; best J drop wins.
#      Every GB greedy steps: NLM Levenberg-Marquardt passes over (position, delta) of the changed set only (FD Jacobian in position,
#      exact in delta; Marquardt diag scaling; projection onto the limits above).
# Efficiency vs v2: X = Ai @ Ac^T cached and rank-updated (rows of changed bodies, columns of replaced candidates).
# Stops at cancel rms <= TGT, stall, or MAXMV (MAXMV = compute budget, 1000 like v2). Checkpoint per cycle (small state; big matrices rebuilt
# on resume). Per-cycle snapshots of the residual map in snaps/. Fails loud on NaN/inf.
import sys,os,pickle,json,argparse,time,numpy as np,healpy as hp
ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=1);ap.add_argument('--NB',type=int,default=3000);ap.add_argument('--NC',type=int,default=20000)
ap.add_argument('--sig',type=float,default=0.078);ap.add_argument('--LC',type=int,default=130);ap.add_argument('--TGTF',default='target_Tobs_l130.npy')
ap.add_argument('--RMAXF',type=float,default=0.975);ap.add_argument('--RMIN',type=float,default=0.01);ap.add_argument('--TOL',type=float,default=1e-4);ap.add_argument('--MAXMV',type=int,default=1000)
ap.add_argument('--GB',type=int,default=25);ap.add_argument('--NLM',type=int,default=3);ap.add_argument('--HFD',type=float,default=0.02);ap.add_argument('--TGT',type=float,default=1.64e-5);ap.add_argument('--STALL',type=float,default=1e-3)
ap.add_argument('--OUT',default='');a=ap.parse_args();O='/workspace/coepdf/joint_solve/phase5/';W=O+'web/inverse/v3/'+a.OUT;os.makedirs(W+'snaps',exist_ok=True)
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit(f'NONFINITE in {n}')
    return x
t00=time.time()
MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;NS=64;T0=2.7255e6;R=21.7439;EPS=0.025*R;LMX=3*NS-1;RHOC=2.775e11*(68.64/100)**2
B=fin(np.load(O+'bg_RA_0.01_St_sm_ns64.npy'),'B');n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T;ll,mm=hp.Alm.getlm(LMX)
Tobs=fin(np.load(O+'web/inverse/v3/'+a.TGTF),'Tobs')
rng=np.random.default_rng(a.seed)
def pos(N):u=rng.normal(size=(N,3));u/=np.linalg.norm(u,axis=1)[:,None];r=rng.uniform((a.RMIN*R)**3,(a.RMAXF*R)**3,N)**(1/3);return u*r[:,None]
y0=pos(a.NB);V=4/3*np.pi*((a.RMAXF*R*1e3)**3-(a.RMIN*R*1e3)**3)/a.NB;MU=RHOC*V   # same draw order as v2
g=rng.normal(size=a.NB);dl0=np.exp(a.sig*g-a.sig**2/2)-1;Yc=pos(a.NC)
def kern(v):
    d=v[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@v/np.linalg.norm(v)**3;return (k-k.mean())*MSUN/KG*RV**2
sel=(ll>=1)&(ll<=a.LC);w=np.where(mm[sel]==0,1.,np.sqrt(2.))
def vec(m):x=hp.map2alm(m,lmax=LMX,iter=1)[sel]*w;return np.r_[x.real,x.imag]
NV=2*int(sel.sum());Ab=vec(B-Tobs)
def proj(v):
    q=np.linalg.norm(v);lo,hi=a.RMIN*R,a.RMAXF*R;return v*np.clip(q,lo,hi)/q
CK=W+f'inv3_s{a.seed}_ckpt.pkl'
y=y0.copy();dl=dl0.copy();moved=np.zeros(a.NB,bool);mch=np.zeros(a.NB,bool);nmv=0;J=[];mu=1e-3;cyc=[];resumed=False
if os.path.exists(CK):
    st=pickle.load(open(CK,'rb'));y,dl,moved,mch,nmv,J,Yc,mu,cyc=[st[k] for k in 'y dl moved mch nmv J Yc mu cyc'.split()];rng.bit_generator.state=st['rng'];resumed=True;print('RESUMED at',nmv,flush=True)
Mb=dl*MU
Ai=np.empty((a.NB,NV));Ac=np.empty((a.NC,NV))
for i in range(a.NB):Ai[i]=vec(kern(y[i]))
print('bodies done',round(time.time()-t00),flush=True)
for c in range(a.NC):Ac[c]=vec(kern(Yc[c]))
print('cands done',round(time.time()-t00),flush=True);fin(Ai,'Ai');fin(Ac,'Ac')
cc=np.einsum('ij,ij->i',Ac,Ac);X=np.empty((a.NB,a.NC))
for b0 in range(0,a.NB,250):X[b0:b0+250]=Ai[b0:b0+250]@Ac.T
fin(X,'X');print('X done',round(time.time()-t00),flush=True)
r=Ab+Mb@Ai;Rmap=B.copy()
for i in range(a.NB):Rmap+=Mb[i]*kern(y[i])
fin(Rmap,'Rmap')
if not resumed:
    J=[float(r@r)];R0=Rmap.copy();np.save(W+f'snaps/s{a.seed}_R0.npy',R0)
    np.savez(W+f'snaps/s{a.seed}_m0000.npz',Rmap=Rmap,y=y,dl=dl,moved=moved,mch=mch,nmv=0)
else:R0=np.load(W+f'snaps/s{a.seed}_R0.npy')
J0B=float(Ab@Ab);print('J_bg',J0B,'J_now',float(r@r),flush=True)
def crms():return float((Rmap-Tobs).std())
def setbody(i,v,m):   # set position v and mass m of body i; exact update of r, Rmap, Ai
    global r,Rmap
    anew=vec(kern(v));r=r+m*anew-Mb[i]*Ai[i];Rmap=Rmap+m*kern(v)-Mb[i]*kern(y[i]);Ai[i]=anew;y[i]=v;Mb[i]=m;dl[i]=m/MU
MMIN=-MU   # delta=-1
while nmv<a.MAXMV:
  ii=np.einsum('ij,ij->i',Ai,Ai)
  for _ in range(a.GB):
    rc=Ac@r;ri=Ai@r;best=(0.,None,None,None)
    base=-2*Mb*ri+Mb**2*ii   # |r - M_i a_i|^2 - |r|^2
    # (b) mass-only
    dm=np.maximum(Mb-ri/ii,MMIN)-Mb;dJb=2*dm*ri+dm**2*ii;kb=int(np.argmin(dJb))
    if dJb[kb]<best[0]:best=(dJb[kb],kb,None,Mb[kb]+dm[kb])
    for b0 in range(0,a.NB,500):   # (a) move + optimal mass
        bs=slice(b0,b0+500);q=rc[None,:]-Mb[bs,None]*X[bs];mp=np.maximum(-q/cc[None,:],MMIN)
        dJ=base[bs,None]+2*mp*q+mp**2*cc[None,:];k=np.unravel_index(np.argmin(dJ),dJ.shape)
        if dJ[k]<best[0]:best=(dJ[k],b0+k[0],k[1],mp[k])
    if best[1] is None or -best[0]<a.TOL*J[-1]:break
    _,i,c,m=best
    if c is None:
        setbody(i,y[i].copy(),m);mch[i]=True
    else:
        mchg=bool(m!=Mb[i]);setbody(i,Yc[c].copy(),m);moved[i]=True;mch[i]|=mchg;Yc[c]=pos(1)[0];Ac[c]=vec(kern(Yc[c]));cc[c]=float(Ac[c]@Ac[c]);X[:,c]=Ai@Ac[c]
    X[i]=Ac@Ai[i];ii[i]=Ai[i]@Ai[i];nmv+=1;J.append(float(r@r));fin(r,'r')
  S=np.where(moved|mch)[0];M0=MU   # LM over (x,y,z,delta) of the changed set
  for _ in range(a.NLM):
    if len(S)==0:break
    Jm=np.empty((NV,4*len(S)))
    for j,i in enumerate(S):
        for k in range(3):e=np.zeros(3);e[k]=a.HFD;Jm[:,4*j+k]=Mb[i]*(vec(kern(y[i]+e))-Ai[i])/a.HFD
        Jm[:,4*j+3]=MU*Ai[i]
    fin(Jm,'Jm');G=Jm.T@Jm;gv=Jm.T@r;dg=np.diag(G).copy();dg=np.maximum(dg,1e-12*dg.max());ok=False;del Jm
    for _t in range(8):
        d=-np.linalg.solve(G+mu*np.diag(dg),gv);old=(r.copy(),Rmap.copy(),Ai[S].copy(),y[S].copy(),Mb[S].copy(),dl[S].copy())
        for j,i in enumerate(S):setbody(i,proj(y[i]+d[4*j:4*j+3]),max(Mb[i]+MU*d[4*j+3],MMIN))
        if float(r@r)<J[-1]:ok=True;mu=max(mu/3,1e-9);break
        r,Rmap=old[0],old[1];Ai[S]=old[2];y[S]=old[3];Mb[S]=old[4];dl[S]=old[5];mu*=10
    J.append(float(r@r));fin(r,'r')
    if not ok:break
  for b0 in range(0,len(S),250):X[S[b0:b0+250]]=Ai[S[b0:b0+250]]@Ac.T
  cr=crms();nd=float(hp.remove_dipole(Rmap-Tobs).std())
  cyc.append(dict(n_moves=nmv,unique_changed=int((moved|mch).sum()),moved=int(moved.sum()),mass_only=int((mch&~moved).sum()),J=J[-1],cancel_rms=cr,resid_nodip_rms=nd,mu=mu,t=round(time.time()-t00)));print('cycle',cyc[-1],flush=True)
  np.savez(W+f'snaps/s{a.seed}_m{nmv:04d}.npz',Rmap=Rmap,y=y,dl=dl,moved=moved,mch=mch,nmv=nmv)
  pickle.dump(dict(y=y,dl=dl,moved=moved,mch=mch,nmv=nmv,J=J,Yc=Yc,mu=mu,cyc=cyc,rng=rng.bit_generator.state),open(CK+'.tmp','wb'));os.replace(CK+'.tmp',CK)
  if cr<=a.TGT:print('TARGET reached',flush=True);break
  if len(cyc)>3 and cyc[-4]['J']-J[-1]<a.STALL*cyc[-4]['J']:print('STALL',flush=True);break
  if len(cyc)>1 and cyc[-1]['n_moves']==cyc[-2]['n_moves'] and cyc[-1]['J']>=cyc[-2]['J']:print('NO PROGRESS',flush=True);break
ch=moved|mch;Mb=dl*MU
Resid=B.copy()
for i in range(a.NB):Resid+=Mb[i]*kern(y[i])
fin(Resid,'Resid');drift=float(np.abs(Resid-Rmap).max());print('incremental-map drift',drift)
ry0=np.linalg.norm(y0,axis=1);ry=np.linalg.norm(y,axis=1)
pin_out=int((ry[ch]>=0.999*a.RMAXF*R).sum());pin_in=int((ry[ch]<=1.001*a.RMIN*R).sum());pin_void=int((dl[ch]<=-1+1e-9).sum())
def Dl(m,L=LMX):l=np.arange(L+1);return l*(l+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(m),lmax=L)*T0**2
def pct(x):return dict(mean=float(x.mean()),median=float(np.median(x)),p10=float(np.percentile(x,10)),p90=float(np.percentile(x,90)),min=float(x.min()),max=float(x.max())) if len(x) else None
dsp=np.linalg.norm(y-y0,axis=1)[moved];dd=(dl-dl0)[ch]
def band(m,l1,l2):f=np.zeros(LMX+1);f[l1:l2+1]=1;return float(hp.alm2map(hp.almxfl(hp.map2alm(m,lmax=LMX,iter=3),f),NS,lmax=LMX).std())
out=dict(args=vars(a),R_Gpc=R,n_bodies=a.NB,MU_Msun_per_unit_delta=MU,n_clusters_start=int((dl0>0).sum()),n_voids_start=int((dl0<0).sum()),n_clusters_end=int((dl>0).sum()),n_voids_end=int((dl<0).sum()),
 delta0_min=float(dl0.min()),delta0_max=float(dl0.max()),delta_end_min=float(dl.min()),delta_end_max=float(dl.max()),
 n_moves=nmv,n_changed_unique=int(ch.sum()),frac_changed=float(ch.mean()),n_moved_position=int(moved.sum()),n_mass_only=int((mch&~moved).sum()),
 J_bg=J0B,J_start=J[0],J_end=J[-1],disp_Gpc=pct(dsp),ddelta_changed=pct(dd),absddelta_changed=pct(np.abs(dd)),
 sign_flips_cluster_void=int((np.sign(dl)!=np.sign(dl0))[ch].sum()),
 rms_nodip_bg=float(hp.remove_dipole(B).std()),rms_nodip_start_sky=float(hp.remove_dipole(R0).std()),rms_nodip_end_sky=float(hp.remove_dipole(Resid).std()),rms_target=float(Tobs.std()),
 resid_nodip_rms_start=float(hp.remove_dipole(R0-Tobs).std()),resid_nodip_rms_end=float(hp.remove_dipole(Resid-Tobs).std()),
 resid_rms_l2_40_end=band(Resid-Tobs,2,40),resid_rms_l41_130_end=band(Resid-Tobs,41,130),resid_rms_l131_191_end=band(Resid-Tobs,131,LMX),
 cancel_rms_end=float((Resid-Tobs).std()),cancel_max_end=float(np.abs((Resid-Tobs)-(Resid-Tobs).mean()).max()),
 r_over_R_changed=dict(final_mean=float((ry[ch]/R).mean()),orig_mean=float((ry0[ch]/R).mean()),frac_final_gt_0p75=float((ry[ch]/R>0.75).mean()),frac_orig_gt_0p75=float((ry0[ch]/R>0.75).mean()),
   frac_uniform_gt_0p75=float(1-((0.75**3-a.RMIN**3)/(a.RMAXF**3-a.RMIN**3)))),
 pinned=dict(at_RMAXF_0p975R=pin_out,at_RMIN=pin_in,at_delta_minus1=pin_void),cycles=cyc,map_drift=drift,moved_absmass_frac_start=float(np.abs(dl0[ch]).sum()/np.abs(dl0).sum()),
 runtime_s=round(time.time()-t00))
json.dump(out,open(W+f'inv3_s{a.seed}.json','w'),indent=1)
np.savez(W+f'inv3_s{a.seed}.npz',y0=y0,y=y,delta0=dl0,delta=dl,M0=dl0*MU,M=Mb,moved=moved,mch=mch,changed=ch,J=np.array(J),Resid=Resid,R0=R0,B=B,Tobs=Tobs,Dl_resid=Dl(Resid),Dl_start=Dl(R0),Dl_obs=Dl(Tobs))
print('OUT',json.dumps({k:v for k,v in out.items() if k!='cycles'}),flush=True)
