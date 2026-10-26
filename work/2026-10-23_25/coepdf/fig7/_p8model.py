# rev7 helper: load the p8 forward model (joint_solve/phase2/p8_field_consistent.py) with its fitted Limber/exact correction CORR.
import sys,json,os,numpy as np
H=os.path.dirname(os.path.abspath(__file__));O='/workspace/coepdf/joint_solve/phase2/'
sys.path.insert(0,O);cwd=os.getcwd()
import p8_field_consistent as P8
os.chdir(cwd)
P8.CORR[:]=np.array(json.load(open(O+'p8_ckpt.json'))['corr'])
J=json.load(open(O+'p8_field_consistent.json'));XB=np.array(J['best_x']);PB=J['best_p']
M,A=P8.M,P8.A
