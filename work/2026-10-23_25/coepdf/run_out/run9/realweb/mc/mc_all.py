import numpy as np,json
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
l=np.arange(2,31);G=np.logspace(-1.2,0.1,25);out=[]
for k in range(len(S)):
  r={}
  for iR in (3,4,5):
    sk=np.array(S[k][iR]);cs=[chi(Dl(sm+e*sk)) for e in G];i=int(np.argmin(cs));Db=Dl(sm+G[i]*sk)[2:31]
    sl=np.polyfit(l[2:],np.log(Db[2:]),1)[0];ov=Db[:3].max()/Pl[:3,1].max()
    r[RS[iR]]=dict(e=float(G[i]),chi2=float(cs[i]),guess=bool(cs[i]<=100 and abs(sl)<0.05))
  out.append(r)
json.dump(out,open('samples/all_guesses.json','w'))
for R in (4.0,8.0,14.0):
  g=[o[R]['guess'] for o in out];c=[o[R]['chi2'] for o in out];print(R,sum(g),np.median(c),min(c))
print('all3',sum(all(o[R]['guess'] for R in o) for o in out))
