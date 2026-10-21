import numpy as np,json
from multiprocessing import Pool
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];pm=np.median(Pl[2:,1]);G=json.load(open('samples/bracket_guesses.json'))
def f(args):
  k,iR=args;sk=np.array(S[k][iR]);med=lambda e:np.median(Dl(sm+e*sk)[4:31])
  st=[]
  for e in E:
    y=Dl(sm+e*sk)[4:31];st.append((np.median(y),abs(np.polyfit(l[2:],np.log(y),1)[0])<0.05))
  i=[i for i in range(9) if st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0]][0]
  a,b=np.log(E[i+1]),np.log(E[i])
  for _ in range(10):
    m=(a+b)/2
    if med(np.exp(m))>pm:b=m
    else:a=m
  e=np.exp((a+b)/2);return e,Dl(sm+e*sk)
out={}
with Pool() as p:
  for R in [4.0,8.0,14.0]:
    iR=RS.index(R);ks=[k for k,o in enumerate(G) if o[str(R)]]
    r=[]
    for j,x in enumerate(p.imap(f,[(k,iR) for k in ks],chunksize=4)):
      r.append(x)
      if j%100==0:print(R,j,flush=True)
    out[R]=(np.array([x[0] for x in r]),np.array([x[1] for x in r]));print(R,len(ks),flush=True)
np.savez('mc_mean_passing.npz',**{f'E{R}':out[R][0] for R in out},**{f'D{R}':out[R][1] for R in out})
