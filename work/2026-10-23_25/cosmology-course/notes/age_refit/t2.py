import sys,time;sys.argv=['x','exp']
exec(open('unified2.py').read().split("if __name__=='__main__'")[0])
for i in range(2):
    t=time.time();e=model(dict(D=58.2,M=1e53,l0=15.3,L=6.07,H=60.,qc=0.23,p=1.0));print(time.time()-t,{k:v for k,v in e.items() if k!='A'} if e else e)
