import numpy as np,json
exec(open('cmbann.py').read().split("meanT=lambda")[0])
sp=json.load(open('split.json')); qi,qo=sp['q_in'],sp['q_out']
V=lambda a,c,t: r_of(c,t)**3-r_of(a,t)**3; e=2e-3
H=lambda a,c: (np.log(V(a,c,1.0))-np.log(V(a,c,1.0-e)))/(3*e)   # (1/3) dlnV/dt at now, code units
Hin,Hout=H(qi,1.0),H(1.0,qo)
json.dump(dict(H_in=Hin,H_out=Hout,H_in_per_Gyr=Hin/t0,H_out_per_Gyr=Hout/t0,ratio_out_in=Hout/Hin),open('splitA.json','w'),indent=1)
print(Hin,Hout,Hout/Hin)
