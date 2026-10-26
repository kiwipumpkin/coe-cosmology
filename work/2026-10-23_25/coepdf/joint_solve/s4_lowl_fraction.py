# Step 3: the remaining low-l excess is the stochastic field itself (D_2 ~ 3150 uK^2). Masses placed as in s2 can cancel a fraction
# (1-f) of the field's own low-l a_lm on a given sky (D_l -> f^2 D_l). Scan f for (a) l=2 only, (b) all l<=12. One extra parameter each.
# Uses the s3 refit field; seed-42 one-sky realisation as in half_figs_fixed.py. Then verifies with masses on an actual seed-42 a_lm sky.
import sys,os,json,numpy as np,healpy as hp
OUT=os.path.abspath('.');os.chdir('/workspace/redweb2/stochfield');sys.path.insert(0,'.')
import acoustic_fit as A
s3=sys.argv[1];z=np.load(os.path.join(OUT,s3));Df,Dfix=z['Df'],z['Dfix_new'];ll=np.arange(len(Df))
g=np.random.default_rng(42).chisquare(2*np.maximum(ll,1)+1)/(2*np.maximum(ll,1)+1)
def st(fl):
    D=Df*fl**2+Dfix;Dr=Df*g*fl**2+Dfix
    return {k:A.stats(X,Df*fl**2,*r)['chi2_planck'] for k,X in (('mean',D),('seed42',Dr)) for r in [(2,2500)]}|{k+'_2_30':A.stats(X,Df*fl**2,2,30)['chi2_planck'] for k,X in (('mean',D),('seed42',Dr))}|{'r2':A.stats(D,Df,2,2500)['r2']}
fs=np.round(np.linspace(0,1,101),2);out={'s3':s3,'none':st(np.ones(len(Df)))}
for name,LC in (('l2only',2),('l_le12',12)):
    sc=[]
    for f in fs:
        fl=np.ones(len(Df));fl[2:LC+1]=f;sc.append((f,st(fl)))
    b=min(sc,key=lambda t:t[1]['mean']);out[name]=dict(best_f=float(b[0]),best=b[1],scan=[[float(f),s['mean'],s['seed42']] for f,s in sc])
json.dump(out,open(os.path.join(OUT,'s4_'+s3.replace('.npz','.json')),'w'),indent=1)
print(json.dumps({k:(v if k in('none','s3') else {'best_f':v['best_f'],**v['best']}) for k,v in out.items()},indent=1))
