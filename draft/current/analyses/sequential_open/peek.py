import sys,pickle,numpy as np
from common import *
import os; os.chdir(SEQ); import scenario_b as sb; os.chdir(OUT)
st=pickle.load(open(sys.argv[1],"rb")); g=st["g"]
X,M,k,left,w=sb.trajectory(g.order)
print(len(g.voids),"voids",len(M),"objs; now %.1f  min %.1f at %d"%(left[-1],left.min(),left.argmin()+1), "web(rms) %.1f"%(rms(g.w)/FLOOR), where(X))
