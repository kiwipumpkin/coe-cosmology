import numpy as np,json
exec(open('/workspace/redweb2/realweb/analyze.py').read().split('O=')[0])
S=np.load('mc_sky.npy',mmap_mode='r');sm=np.load('/workspace/redweb2/sm_red.npy');RS=[0.5,0.9,2.0,4.0,8.0,14.0]
l=np.arange(2,31);E=[1,.5,.2,.1,.05,.02,.01,.005,.002,.0005];pm=np.median(Pl[2:,1]);out=[]
for k in range(len(S)):
  r={}
  for iR,R in enumerate(RS):
    sk=np.array(S[k][iR]);st=[]
    for e in E:
      y=Dl(sm+e*sk)[4:31];st.append((np.median(y),abs(np.polyfit(l[2:],np.log(y),1)[0])<0.05))
    ok=any(st[i][1] and st[i+1][1] and st[i][0]>=pm>=st[i+1][0] for i in range(len(E)-1))
    r[R]=bool(ok)
  out.append(r)
  if k%100==0:print(k,flush=True)
json.dump(out,open('samples/bracket_guesses.json','w'))
for R in RS:print('R',R,sum(o[R] for o in out))
print('all 4/8/14',sum(o[4.0] and o[8.0] and o[14.0] for o in out))
print('ex_pass',[k for k,o in enumerate(out) if o[4.0] and o[8.0] and o[14.0]][:5],'ex_fail',[k for k,o in enumerate(out) if not(o[4.0] or o[8.0] or o[14.0])][:5])
