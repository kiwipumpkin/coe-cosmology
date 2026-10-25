from core import *
w=json.load(open("../../migkas-geometry-fit.json"))["wls_fit"]
t=time.time(); s=ltb(w["D"],w["mcc"],w["kin"],w["kout"],True); print("t",time.time()-t)
c=shifted(s["af"],w["l0"]); print("chi",chi(c),r2(c),"Ht",s["Ht"],"E+",(s["E"]>0).sum())
t=time.time(); n=nsky(w["D"],w["mcc"],w["kin"],w["kout"]); print("newt t",time.time()-t, chi(shifted(n,w["l0"])))
