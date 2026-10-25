import numpy as np,json
src=open('jk.py').read().split("allk=np.ones")[0];exec(src)
allk=np.ones(len(z),bool);A0,B0=curves(allk);la,Bg=fit(allk);res=y-(la+Bg*x)
out={'grid':grid.tolist(),'A0':A0.tolist(),'B0':B0.tolist(),'AF':AFm.tolist(),'BL':BLm.tolist()};rep={}
rng=np.random.default_rng(2)
for tag,lo,hi,key in (('Y',15,80,'AB'),('P',260,335,'B')):
  S=(grid>=lo)&(grid<=hi);ref=(AFm if key=='A' else BLm)
  dev=lambda A1,B1:((A1-AFm)[S].mean(),(B1-BLm)[S].mean())
  d0=dev(A0,B0);c=(lo+hi)/2;hw=(hi-lo)/2+45
  idx=np.where(np.abs((l-c+180)%360-180)<=hw)[0];rows=[]
  for i in idx:
    k=allk.copy();k[i]=False;A1,B1=curves(k);d=dev(A1,B1)
    if key=='AB':sc=-(d[0]-d0[0])/0.010+(d[1]-d0[1])/0.009
    else:
      j=0 if key=='A' else 1;sc=-(d[j]-d0[j])*np.sign(d0[j])
    rows.append((sc,int(i)))
  rows.sort(key=lambda r:-r[0]);top=[r[1] for r in rows]
  R={'N':len(idx),'dev0':d0,'list10':[(name[i].strip(),round(l[i],1),round(b[i],1),round(z[i],3),T[i],round(res[i],2)) for i in top[:10]]}
  for kk in (3,5,10):
    k=allk.copy();k[top[:kk]]=False;A1,B1=curves(k);R[f'top{kk}']=dev(A1,B1);out[f'{tag}A{kk}']=A1.tolist();out[f'{tag}B{kk}']=B1.tolist()
  v=[]
  for _ in range(100):
    k=allk.copy();k[rng.choice(idx,10,replace=False)]=False;v.append(dev(*curves(k)))
  v=np.array(v);R['rand10_mean']=v.mean(0).tolist();R['rand10_sd']=v.std(0).tolist();rep[tag]=R
json.dump(out,open('curves.json','w'));json.dump(rep,open('jk2.json','w'),default=float,indent=1)
print(json.dumps(rep,default=float,indent=1))
