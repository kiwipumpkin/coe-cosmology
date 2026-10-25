import numpy as np,json,re
src=open('jk.py').read().split("allk=np.ones")[0]
exec(src)
allk=np.ones(len(z),bool);A0,B0=curves(allk);la,Bg=fit(allk);res=y-(la+Bg*x)
rA=A0-AFm;rB=B0-BLm
for nmw,(lo,hi) in (('30-70',(30,70)),('210-250',(210,250)),('40-60',(40,60)),('220-240',(220,240))):
  s=(grid>=lo)&(grid<=hi);print(nmw,'meandev A %+.4f  B %+.4f'%(rA[s].mean(),rB[s].mean()))
S=[np.where(grid==c)[0][0] for c in (220,230,240)];bA0=rA[S].mean();bB0=rB[S].mean();print('220-240 core A %+.4f B %+.4f'%(bA0,bB0))
idx=np.where(np.abs((l-230+180)%360-180)<=65)[0];rows=[]
for i in idx:
  k=allk.copy();k[i]=False;A1,B1=curves(k);dA=(A1-A0)[S].mean();dB=(B1-B0)[S].mean()
  rows.append((-(dA*np.sign(bA0))/abs(bA0)-(dB*np.sign(bB0))/abs(bB0),name[i].strip(),l[i],b[i],z[i],T[i],NH[i],res[i],dA,dB,int(i)))
rows.sort(key=lambda r:-r[0]);print('N window',len(idx))
for k_,r in enumerate(rows[:10]):print(k_+1,r[1],'%.1f %+.1f %.3f %.1f NH%.1f res%+.2f dA%+.4f dB%+.4f'%r[2:10])
for kk in (5,10):
  k=allk.copy();k[[r[10] for r in rows[:kk]]]=False;A1,B1=curves(k);print('top',kk,'A %+.4f B %+.4f'%((A1-AFm)[S].mean(),(B1-BLm)[S].mean()))
