import sys,pickle; sys.path.insert(0,'.'); from runall import *
def job(c):
    to,mu,u,v0=c
    try: w=weak_coast(mu,u,v0,to)
    except Exception as e: w=None
    return c,w,gr_coast(mu,u,v0,to)
if __name__=='__main__':
    T=[float(a) for a in sys.argv[1:]]
    C=[(to,mu,u,v0) for to in T for v0 in (0.1,0.3,0.5) for u in (0,0.9999) for mu in (0.0,0.1)]
    with mp.Pool(8) as P: R=dict((c,(w,g)) for c,w,g in P.map(job,C))
    pickle.dump(R,open('speeds_t_%s.pkl'%'_'.join(sys.argv[1:]),'wb'))
    for c,(w,g) in sorted(R.items()):
        nb=int(np.isnan(g[0]).sum()); s=lambda r:"%.4f %.4f %.4f"%(stats(r)[0],stats(r)[1],stats(r)[3])
        print(c,"weak",s(w) if w is not None else "BROKE","| GR",s(norm(g[0])) if nb==0 else "missing %d"%nb,"hz",g[2],flush=True)
