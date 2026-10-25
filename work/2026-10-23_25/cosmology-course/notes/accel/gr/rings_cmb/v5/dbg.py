import numpy as np,v5x as V
from scipy.integrate import solve_ivp
M=0.008532635720463873;Rrec=0.1
wo,t0=V.setup(M);wmin=-np.sqrt(2*M/Rrec)*0.9995;S,T,W=V.build(M,t0,wmin)
# R along the CoE-ward axis: at each tau, minimum R over labels
for t in [0.9,0.7,0.5,0.4,0.3,0.2,0.1]:
    print(t,[round(float(S(t,w)[0,0]),3) for w in [wmin,-0.3,0,0.3,0.6,1.0]])
