import sys,os,json,numpy as np
sys.argv=['x','exp']; os.environ['RING']='low'
exec(open('unified2.py').read().split("if __name__=='__main__'")[0])
from scipy.stats import chi2 as C2, norm
z=json.load(open('s12_exp_ckpt0.json'))['x']
P=dict(D=float(np.exp(z[0])),M=float(np.exp(z[1])),l0=z[2]%360,L=float(np.exp(z[3])),p=z[4],H=70.0)
out=dict(params=P)
mig=model(dict(P,qc=0.3),True); out['mig']=dict(chi2=mig['chi2'],R2=mig['R2'],A=mig['A'])
d=45.1366-mig['chi2']; out['mig']['sigma_k3']=float(norm.isf(C2.sf(d,3))); out['mig']['sigma_k5']=float(norm.isf(C2.sf(d,5)))
for tag,sr in (('grav',False),('sr_only',True)):
    e=model(dict(P),sr); e.pop('A'); out['H70_'+tag]=e
H13=70*out['H70_grav']['age']/13
for tag,sr in (('grav',False),('sr_only',True)):
    e=model(dict(P,H=H13,M=P['M']*(H13/70)**2),sr); e.pop('A'); e['H']=H13; out['age13_'+tag]=e
json.dump(out,open('step2_best.json','w'),indent=1,default=float)
print(json.dumps({k:v for k,v in out.items() if k!='mig'},default=float,indent=0)[:3000]); print(out['mig']['chi2'],out['mig']['R2'],out['mig']['sigma_k3'],out['mig']['sigma_k5'])
