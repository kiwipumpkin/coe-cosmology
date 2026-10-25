import sys; sys.path.insert(0,'.'); import numpy as np, iso
for m in (0.1,0.3):
    Mt=iso.Metric(m); Xo=np.array([2.0,0]); ems=[]
    for t in np.radians(np.arange(0,360,45)):
        X=Xo+np.array([np.cos(t),np.sin(t)]); ems.append(((lambda X:(lambda s:(X,np.zeros(2))))(X),-6.0))
    r,_=iso.observe(Mt,Xo,np.zeros(2),ems,0.0)
    ex=np.array([np.sqrt(Mt.AB(0,em(0)[0])[0]/Mt.AB(0,Xo)[0]) for em,_ in ems]); print("static iso m",m,"max err %.1e"%np.max(abs(r/ex-1)))
