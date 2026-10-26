# Interim evaluation of the best finished p4 multistart point (from p4_ckpt.json) with the p4 report(); no corr refresh / polish.
import sys,json,numpy as np
sys.argv=[sys.argv[0],'/workspace/coepdf/joint_solve/s3_s2_s2000_N3000_L12_M3e+16_r3.5_lam1e-05.npz']
exec(open('/workspace/coepdf/joint_solve/phase2/p4_tt_extended.py').read().split("if __name__")[0])
R=json.load(open(OUT+'p4_ckpt.json'));R.sort(key=lambda t:t[2]);x0=np.r_[x12,0,0,2.0]
b=report(x0)[0];rep,D,Dr,Df=report(np.array(R[0][1]))
out=dict(n_starts_done=len(R),best_fit_measure=R[0][2],baseline_fit_measure=b['fit_measure'],best_x=R[0][1],best_p=P(np.array(R[0][1])),best=rep,baseline=b)
json.dump(out,open(OUT+'p4b_interim.json','w'),indent=1,default=float);np.savez(OUT+'p4b_interim.npz',D=D,Dr=Dr,Df=Df,Dfix=Dfix)
for k in('mean','seed42'):print(k,'base',{r:round(b[k][r]['chi2_planck'],2) for r in b[k]},'best',{r:round(rep[k][r]['chi2_planck'],2) for r in rep[k]})
print('fit',b['fit_measure'],R[0][2],len(R))
