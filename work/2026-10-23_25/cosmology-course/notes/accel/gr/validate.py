from grcore import *
# 1) static emitters on the ring, static observer: ratio must equal sqrt(g_tt(e)/g_tt(o))
for M in (0.05,0.1,0.3):
    em=[]
    for t in th[::30]:
        X=np.array([np.cos(t),np.sin(t)]); em.append((lambda X:(lambda s:(X,fourvel(M,*polar(X),np.zeros(2)))))(X))
    Xo=np.zeros(2); Uo=fourvel(M,*polar(Xo),np.zeros(2))
    rat,bad=observe(M,Xo,Uo,em,0.0,-6.0)
    exact=np.array([np.sqrt((1-2*M/polar(np.array([np.cos(t),np.sin(t)]))[0])/(1-2*M/2)) for t in th[::30]])
    print("static M",M,"max |GR/exact-1| = %.2e"%np.max(np.abs(rat/exact-1)),"bad",bad)
# 2) null condition check & flat-space Doppler check: moving emitter, M=0
