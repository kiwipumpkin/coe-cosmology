import sys,json,shutil;sys.path.insert(0,'.');from common import *
J='/workspace/coepdf/joint_solve/phase5/planarity/planarity.json';shutil.copy(J,J+'.bak');R=json.load(open(J));b=R['CoE_bg_only'];ax={}
for l in range(2,LM+1):
    ev,e=eig(ALM['CoE_bg_only'],l);n=e[:,0];ax[l]=n;x=lb(n)
    b[str(l)]['best_axis_lb_OLD_degenerate']=b[str(l)]['best_axis_lb'];b[str(l)]['ang_to_AoE_OLD']=b[str(l)]['ang_to_AoE']
    b[str(l)].update(best_axis_lb=list(x),ang_to_AoE=ang(n,naoe),axis_type='min-eigenvector (symmetry axis); field axisymmetric',eig_norm=ev.tolist(),ang_to_CoE=ang(n,nco))
b['ang_l2_l3_OLD']=b['ang_l2_l3'];b['ang_l2_l3']=ang(ax[2],ax[3])
b['FLAG']='2026-10-27: axis entries recomputed as min-eigenvector of <L_i L_j>; old max-eigenvector axes were degenerate for an m=0 field. p_Lz2/p_ang values for this sky are not meaningful.'
json.dump(R,open(J,'w'),indent=1);print(json.dumps(b,indent=1))
