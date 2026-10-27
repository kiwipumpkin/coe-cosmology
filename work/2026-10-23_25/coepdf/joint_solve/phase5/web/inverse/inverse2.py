# Inverse problem (Chris, Oct 26): start from a statistically UNIFORM random 3-D web of clusters+voids (uniform in volume out to the
# CMB shell, physical lognormal densities so delta>=-1 by construction, zero mean delta), on the Route A eps=0.01 St background.
# Solve for the MINIMAL set of bodies that must be displaced (masses fixed) so bg+web reproduces the Planck low-l sky (same target
# Tobs and weighted alm objective as p5_comp5.py: l=1..LC fitted, l=LC+1..LP penalised with weight WH).
# v2 (inverse2.py): greedy batches of GB moves alternated with NLM Levenberg-Marquardt refinements of the positions of the moved set only
# (finite-difference Jacobian, h=HFD Gpc); stops when cancel rms |Resid-Tobs| <= TGT (default 2.55e-5 = old free-mass web level) or stall.
# Method: greedy sparse moves (L0 on count): each step moves the single body to the single candidate position (from NC random
# candidates uniform in the same volume) that most lowers J; exact rank-1 update. Stops when relative J drop < TOL per move, or MAXMV.
# Only geometric limits: r in [RMIN,RMAXF]*R (RMAXF=0.975 = one softening length inside shell; RMIN=0.01 avoids r=0). Fails on NaN.
import sys,os,pickle,json,argparse,numpy as np,healpy as hp
ap=argparse.ArgumentParser();ap.add_argument('--seed',type=int,default=1);ap.add_argument('--NB',type=int,default=3000);ap.add_argument('--NC',type=int,default=20000)
ap.add_argument('--sig',type=float,default=0.078);ap.add_argument('--LC',type=int,default=40);ap.add_argument('--LP',type=int,default=60);ap.add_argument('--WH',type=float,default=1.0)
ap.add_argument('--RMAXF',type=float,default=0.975);ap.add_argument('--RMIN',type=float,default=0.01);ap.add_argument('--TOL',type=float,default=1e-4);ap.add_argument('--MAXMV',type=int,default=20000)
ap.add_argument('--GB',type=int,default=25);ap.add_argument('--NLM',type=int,default=3);ap.add_argument('--HFD',type=float,default=0.02);ap.add_argument('--TGT',type=float,default=2.55e-5);ap.add_argument('--STALL',type=float,default=1e-3)
a=ap.parse_args();O='/workspace/coepdf/joint_solve/phase5/';W=O+'web/inverse/'
def fin(x,n):
    if not np.all(np.isfinite(x)):raise SystemExit(f'NONFINITE in {n}')
    return x
MSUN=1.989e30;KG=7.7e51/0.0039;RV=14.3;NS=64;T0=2.7255e6;R=21.7439;EPS=0.025*R;LMX=3*NS-1;RHOC=2.775e11*(68.64/100)**2
B=fin(np.load(O+'bg_RA_0.01_St_sm_ns64.npy'),'B');n=np.array(hp.pix2vec(NS,np.arange(hp.nside2npix(NS)))).T;ll,mm=hp.Alm.getlm(LMX)
lw=np.arange(LMX+1).astype(float);S_=hp.ud_grade(hp.read_map('/workspace/planck/smica.fits'),NS)/2.7255;as_=hp.map2alm(S_,lmax=LMX,iter=3)
wt=np.where(lw<=20,1.,np.where(lw<40,np.cos(0.5*np.pi*(lw-20)/20)**2,0.));wt[:2]=0;Tobs=hp.alm2map(hp.almxfl(as_,wt),NS,lmax=LMX)
rng=np.random.default_rng(a.seed)
def pos(N):u=rng.normal(size=(N,3));u/=np.linalg.norm(u,axis=1)[:,None];r=rng.uniform((a.RMIN*R)**3,(a.RMAXF*R)**3,N)**(1/3);return u*r[:,None]
y0=pos(a.NB);V=4/3*np.pi*((a.RMAXF*R*1e3)**3-(a.RMIN*R*1e3)**3)/a.NB
g=rng.normal(size=a.NB);dl=np.exp(a.sig*g-a.sig**2/2)-1;# lognormal: delta>-1, E[delta]=0; sig default 0.078 = std of delta the old free-mass web needed
Mb=dl*RHOC*V;Yc=pos(a.NC)
def kern(v):
    d=v[None,:]-R*n;dd=(np.sum(d*d,1)+EPS**2)**1.5;k=np.sum(n*d,1)/dd-n@v/np.linalg.norm(v)**3;return (k-k.mean())*MSUN/KG*RV**2
sel=((ll>=1)&(ll<=a.LC))|((ll>a.LC)&(ll<=a.LP));w=np.where(mm[sel]==0,1.,np.sqrt(2.))*np.where(ll[sel]>a.LC,a.WH,1.)
def vec(m):x=hp.map2alm(m,lmax=LMX,iter=1)[sel]*w;return np.r_[x.real,x.imag]
Ab=vec(B-Tobs);Ai=np.array([vec(kern(v)) for v in y0],dtype=np.float64);print('bodies done',flush=True)
Ac=np.array([vec(kern(v)) for v in Yc],dtype=np.float64);print('cands done',flush=True);fin(Ai,'Ai');fin(Ac,'Ac')
cc=np.einsum('ij,ij->i',Ac,Ac).astype(np.float64);r=Ab+Mb@Ai;J0B=float(Ab@Ab);J=[float(r@r)];y=y0.copy();moved=np.zeros(a.NB,bool);nmv=0;hist=[]
print('J_bg',J0B,'J_start(bg+random web)',J[0],flush=True)
Rmap=B+np.array([kern(v) for v in y0]).T@Mb;fin(Rmap,'Rmap');R0=Rmap.copy();mu=1e-3;cyc=[];npin=0
def crms():return float((Rmap-Tobs).std())
def setpos(i,v):
    global r,Rmap
    anew=vec(kern(v));r=r+Mb[i]*(anew-Ai[i]);Rmap=Rmap+Mb[i]*(kern(v)-kern(y[i]));Ai[i]=anew;y[i]=v
def proj(v):
    q=np.linalg.norm(v);lo,hi=a.RMIN*R,a.RMAXF*R;return (v*np.clip(q,lo,hi)/q,bool(q<lo or q>hi))
Jc0=J[-1];CK=W+f'inv_s{a.seed}_ckpt.pkl'
if os.path.exists(CK):
    st=pickle.load(open(CK,'rb'));y,moved,nmv,J,Yc,Ac,cc,Ai,r,Rmap,mu,cyc=[st[k] for k in 'y moved nmv J Yc Ac cc Ai r Rmap mu cyc'.split()];rng.bit_generator.state=st['rng'];print('RESUMED at',nmv,flush=True)
while nmv<a.MAXMV:
  ii=np.einsum('ij,ij->i',Ai,Ai)
  for _ in range(a.GB):
    rc=Ac@r;ri=Ai@r;best=(0,None,None)
    for b0 in range(0,a.NB,500):
        bs=slice(b0,b0+500);X=Ai[bs]@Ac.T;m=Mb[bs][:,None]
        dJ=2*m*(rc[None,:]-ri[bs,None])+m**2*(cc[None,:]-2*X+ii[bs,None]);k=np.unravel_index(np.argmin(dJ),dJ.shape)
        if dJ[k]<best[0]:best=(dJ[k],b0+k[0],k[1])
    if best[1] is None or -best[0]<a.TOL*J[-1]:break
    i,c=best[1],best[2];setpos(i,Yc[c].copy());ii[i]=Ai[i]@Ai[i];moved[i]=True;nmv+=1;J.append(float(r@r));fin(r,'r')
    Yc[c]=pos(1)[0];Ac[c]=vec(kern(Yc[c]));cc[c]=float(Ac[c]@Ac[c])
  S=np.where(moved)[0]
  for _ in range(a.NLM):   # LM on moved positions
    Jm=np.empty((len(r),3*len(S)))
    for j,i in enumerate(S):
        for k in range(3):e=np.zeros(3);e[k]=a.HFD;Jm[:,3*j+k]=Mb[i]*(vec(kern(y[i]+e))-Ai[i])/a.HFD
    G=Jm.T@Jm;gv=Jm.T@r;ok=False
    for _t in range(8):
        d=-np.linalg.solve(G+mu*np.trace(G)/len(G)*np.eye(len(G)),gv);old=(r.copy(),Rmap.copy(),Ai[S].copy(),y[S].copy())
        for j,i in enumerate(S):setpos(i,proj(y[i]+d[3*j:3*j+3])[0])
        if float(r@r)<J[-1]:ok=True;mu=max(mu/3,1e-9);break
        r,Rmap=old[0],old[1];Ai[S]=old[2];y[S]=old[3];mu*=10
    J.append(float(r@r));fin(r,'r')
  cr=crms();cyc.append(dict(n_moves=nmv,unique=int(moved.sum()),J=J[-1],cancel_rms=cr,mu=mu));print('cycle',cyc[-1],flush=True)
  pickle.dump(dict(y=y,moved=moved,nmv=nmv,J=J,Yc=Yc,Ac=Ac,cc=cc,Ai=Ai,r=r,Rmap=Rmap,mu=mu,cyc=cyc,rng=rng.bit_generator.state),open(CK+'.tmp','wb'));os.replace(CK+'.tmp',CK)
  if cr<=a.TGT:print('TARGET reached',flush=True);break
  if len(cyc)>3 and cyc[-4]['J']-J[-1]<a.STALL*cyc[-4]['J']:print('STALL',flush=True);break
npin=int(sum(proj(y[i])[1] or np.linalg.norm(y[i])>=0.999*a.RMAXF*R for i in np.where(moved)[0]))
Resid=B+np.array([kern(v) for v in y]).T@Mb;fin(Resid,'Resid');drift=float(np.abs(Resid-Rmap).max());print('incremental-map drift',drift)
def nodip(m):return float(hp.remove_dipole(m).std())
def Dl(m,L=LMX):l=np.arange(L+1);return l*(l+1)/(2*np.pi)*hp.anafast(hp.remove_dipole(m),lmax=L)*T0**2
dsp=np.linalg.norm(y-y0,axis=1)[moved];ry0=np.linalg.norm(y0,axis=1);ry=np.linalg.norm(y,axis=1)
out=dict(args=vars(a),R_Gpc=R,n_bodies=a.NB,n_clusters=int((Mb>0).sum()),n_voids=int((Mb<0).sum()),delta_min=float(dl.min()),delta_max=float(dl.max()),
 n_moves=nmv,n_moved_unique=int(moved.sum()),frac_moved=float(moved.mean()),J_bg=J0B,J_start=J[0],J_end=J[-1],
 disp_Gpc=dict(mean=float(dsp.mean()),median=float(np.median(dsp)),p10=float(np.percentile(dsp,10)),p90=float(np.percentile(dsp,90)),max=float(dsp.max())) if len(dsp) else None,
 rms_nodip_bg=nodip(B),rms_nodip_start=nodip(R0),rms_nodip_end=nodip(Resid),rms_nodip_Tobs=nodip(Tobs),cancel_rms_end=float((Resid-Tobs).std()),cancel_max_end=float(np.abs((Resid-Tobs)-(Resid-Tobs).mean()).max()),
 cancel_rms_start=float((R0-Tobs).std()),moved_r_over_R=dict(final_mean=float((ry[moved]/R).mean()),orig_mean=float((ry0[moved]/R).mean()),frac_final_gt_0p75=float((ry[moved]/R>0.75).mean()),frac_uniform_gt_0p75=float(1-((0.75**3-a.RMIN**3)/(a.RMAXF**3-a.RMIN**3)))),
 n_pinned_at_RMAXF_or_RMIN=npin,cycles=cyc,map_drift=drift,frac_final_within_1pct_of_RMAXF=float((ry[moved]>0.99*a.RMAXF*R).mean()),moved_mass_frac=float(np.abs(Mb[moved]).sum()/np.abs(Mb).sum()))
json.dump(out,open(W+f'inv_s{a.seed}.json','w'),indent=1)
np.savez(W+f'inv_s{a.seed}.npz',y0=y0,y=y,M=Mb,delta=dl,moved=moved,J=np.array(J),Resid=Resid,R0=R0,B=B,Tobs=Tobs,Dl_resid=Dl(Resid),Dl_start=Dl(R0),Dl_obs=Dl(Tobs))
print('OUT',json.dumps(out),flush=True)
