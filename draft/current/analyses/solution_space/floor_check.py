import json,numpy as np,sys
src=open("analyze.py").read().split("res={}")[0]; exec(src)
ext=json.load(open("extend_mixed.json")); fl=ext[-1]
sols=[dict(fl,src="warm")]
for f in ("cold_CCCCCCVVVVV_polished.json","cold_CCCCCCVVVVVV_polished.json"):
    for r in json.load(open(f)): sols.append(dict(r,src=f.split("_")[1]))
for thr in (1.0,2.0):
    good=sorted([s for s in sols if s["left"]<thr],key=lambda s:s["left"]); reps=[]
    for s in good:
        c=canon(s,s["kinds"]); c["src"]=s["src"]
        for rp in reps:
            if same(rp,c): rp["n"]+=1; break
        else: c["n"]=1; reps.append(c)
    print(f"<{thr}x: solutions {len(good)}, distinct {len(reps)}:",[(round(r['left'],2),r['src'],r['n'],len(r['k'])) for r in reps])
k=np.array(list(fl["kinds"])); X=np.array(fl["X"]); M=np.array(fl["M"]); Rv=np.array(fl["Rv"]); d=np.linalg.norm(X,axis=1)
print("floor",fl["kinds"],round(fl["left"],3),"amps",np.round(fl["amps"],2))
for i in range(len(k)): print(k[i],"d=%.1f l=%.0f b=%.0f M=%.2e Rv=%.2f dCoE=%.1f"%(d[i],np.degrees(np.arctan2(X[i,1],X[i,0]))%360,np.degrees(np.arcsin(X[i,2]/d[i])),M[i],Rv[i],np.linalg.norm(X[i]-C_G)))
print("clusters %.2e voids %.2e net %.2e"%(M[k=="C"].sum(),M[k=="V"].sum(),M.sum()))
for s in ext: print(s["N"],round(s["left"],2),"top m",(np.argsort(s["amps"])[::-1][:3]+1).tolist(),"caps",s.get("caps_mass"),s.get("caps_void"))
