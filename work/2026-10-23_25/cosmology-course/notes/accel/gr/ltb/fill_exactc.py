import pickle,sys; sys.argv=['x']; from run_exactc import *
def runfill(core,v0,T,old):
    L=LTB(core_frac=core,rmax=150.0,nr=900); g=1/np.sqrt(1-v0*v0); ems=[]; l0=[]
    for t in th:
        e=np.array([np.cos(t),np.sin(t)]); ts=L.t1+g*v0*(1+e[0]); w=relsum(v0,1.0,e); w/=np.linalg.norm(w)
        at,_=worldline(L,C0+np.array([g*(e[0]+v0*v0),e[1]]),w,ts,T); ems.append((at,ts)); l0.append(g*(1+v0*e[0]))
    oat,_=worldline(L,C0+np.array([g*v0*v0,0]),np.array([v0,0]),L.t1+g*v0,T); y=oat(T); K=L.kdu
    out=old.copy()
    for i in np.where(np.isnan(old))[0]:
        for j in [k for k in (i-1,i+1,i-2,i+2) if 0<=k<len(th) and not np.isnan(old[k])]:
            pair=[ems[j],ems[i]]; cand=[j,i]
            def fk(t,r,vph,lnk,vp):
                if t==T and lnk==0.0: return K(t,r,vph,lnk,vp)
                q=min(cand,key=lambda c: abs(ems[c][0](t)[0]-r) if ems[c][0](t) is not None else 1e9); e=ems[q][0](t)
                return kdot(L,t,r,vph,lnk,e[2:4],e[4]+np.log(l0[q]))
            L.kdu=fk; r=observe(L,y,pair,T); L.kdu=K
            if not np.isnan(r[1]): out[i]=r[1]; print('filled',core,v0,i,'check neighbour',r[0],old[j]); break
    return out
d=pickle.load(open('exactc.pkl','rb'))
for a,x in list(d.items()):
    if np.isnan(x).any(): d[a]=runfill(a[0],a[1],2.2222222222222223+a[2],x); print(a,np.isnan(d[a]).sum(),np.array(stats(norm(d[a])))[[0,1,3]])
pickle.dump(d,open('exactc.pkl','wb'))
