import numpy as np
exec(open('maxwin.py').read().split('# --- family A')[0])
mus=np.concatenate([np.linspace(0,2,11),np.linspace(2.5,12,20)]); taus=np.linspace(0,0.9,19); res={}
for RD in [0.05,0.246]:
  W=np.full((len(taus),len(mus)),np.nan); V=W.copy()
  for i,ta in enumerate(taus):
    for j,mu in enumerate(mus):
      o=evalf(*fieldA(mu,ta),RD)
      if o and o["ok"]: W[i,j]=o["max"]["win"]; V[i,j]=o["vmax"]
  res[RD]=(W,V); i,j=np.unravel_index(np.nanargmax(W),W.shape); print(RD,taus[i],mus[j],W[i,j],V[i,j],flush=True)
  Wc=np.where(V<1,W,np.nan); i,j=np.unravel_index(np.nanargmax(Wc),W.shape); print(' v<c',taus[i],mus[j],Wc[i,j],flush=True)
np.save("maxwin_A2.npy",dict(mus=mus,taus=taus,res=res),allow_pickle=True)
