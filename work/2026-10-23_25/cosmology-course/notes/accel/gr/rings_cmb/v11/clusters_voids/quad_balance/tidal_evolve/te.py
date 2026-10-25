import sys,os,json,numpy as np
from scipy.spatial import cKDTree
sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir('/workspace/clump_scan/sequential');import seq_core as S;os.chdir(H)
th=np.radians(np.arange(360)) if not hasattr(S,'TH') else S.TH;F=S.FLOOR;G=6.674e-11
def mm(v):m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/F
f=np.load('/workspace/clump_scan/sequential_open/final_comp_s3.npz');X0=f['Xp'].astype(float);r=f['rp'].astype(float);M=f['Mp'].astype(float)
# CoE: point mass at C_G whose pattern reproduces the smooth ring R (fit the one amplitude)
MC=float(os.environ.get('MCOE','1e54'));fitq=float('nan')  # CoE point mass at C_G; only matters for both_raw
print('CoE mass fitted to R: %.3e kg, rel. misfit %.2e'%(MC,fitq))
GPC=S.GPC
def g_coe(X):   # peculiar (tidal) acceleration relative to us: g(x)-g(0), m/s^2
    d=S.C_G-X;d0=S.C_G;return G*MC*(d/np.linalg.norm(d,axis=1)[:,None]**3-d0/np.linalg.norm(d0)**3)/GPC**2
def g_self(X):  # web self-gravity, Plummer-softened by body radius (bodies are extended), signed masses (voids repel)
    out=np.zeros_like(X)
    for i in range(0,len(X),500):
        d=X[None,:,:]-X[i:i+500,None,:];e2=(r[None,:]**2+r[i:i+500,None]**2)
        q=(d*d).sum(-1)+e2;out[i:i+500]=G*np.einsum('ijk,j->ik',d/q[...,None]**1.5,M)/GPC**2
    return out
def gu(X):return G*np.einsum('ijk,j->ik',(0-X)[:,None,:]*0+(X[None,:,:]-0)*0,M) # unused
gc=g_coe(X0);gs=g_self(X0);gs-=0  # (our own pull on us cancels in ring kernel's g(0) term; displacements of bodies only)
def resid(X):return S.R-S.pattern(X,M)
base=resid(X0);a0=mm(base);out=dict(M_coe=MC,start=dict(total=float(S.rms(base)/F),per_m=a0[:12].tolist()))
print('start total %.2f m2 %.2f'%(out['start']['total'],a0[1]))
tree=cKDTree(X0)
def crossings(X):  # overlapping body pairs (shells start to cross)
    t=cKDTree(X);p=t.query_pairs(2*r.max());p=np.array(list(p))
    if len(p)==0:return 0
    return int(np.sum(np.linalg.norm(X[p[:,0]]-X[p[:,1]],axis=1)<r[p[:,0]]+r[p[:,1]]))
n0=crossings(X0);print('overlaps at start',n0)
nt=lambda g:g/np.sqrt(np.mean(np.sum(g*g,1)))
for name,g in (('tidal',gc),('self',gs),('both_raw',gc+gs),('mix25',0.25*nt(gc)+0.75*nt(gs)),('mix50',0.5*nt(gc)+0.5*nt(gs)),('mix75',0.75*nt(gc)+0.25*nt(gs))):
    gr=np.sqrt(np.mean(np.sum(g*g,1)));out[name]=[]
    for s in (0.01,0.03,0.1,0.2,0.3,0.5,1.0,2.0):   # rms displacement, Gpc  (D = s/rms|g|)
        X=X0+s/gr*g;inb=(np.linalg.norm(X,axis=1)<=30)&(np.linalg.norm(X-S.C_G,axis=1)>=S.A_IN_G)
        v=resid(X);a=mm(v);row=dict(rms_disp_Gpc=s,D=s/gr,total=float(S.rms(v)/F),per_m=a[:12].tolist(),m2=float(a[1]),m13_60=float(np.sqrt(np.sum(a[12:]**2))),overlaps=crossings(X)-n0,out_of_region=int((~inb).sum()))
        out[name].append(row);print(name,s,'total %.2f m2 %.2f m1 %.1f m3-12max %.1f m13+ %.1f new overlaps %d outside %d'%(row['total'],a[1],a[0],np.max(a[2:12]),row['m13_60'],row['overlaps'],row['out_of_region']),flush=True)
json.dump(out,open('te.json','w'),indent=1)
