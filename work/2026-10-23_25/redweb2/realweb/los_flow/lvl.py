exec(open('an.py').read().split('SR=json')[0])
import json;out={}
for f in sorted(glob.glob(O+'comp/u*.npz')):
    k=f.split('/u')[1][:-4];z=dict(np.load(f));s=dict(np.load(ST+f'comp/u{k}.npz'));r={}
    for v in ('grow','fixed'):
        D=[Dl(sky(z,s,e,v)) for e in E];md=np.array([np.median(x[4:31]) for x in D]);sl=[np.polyfit(l[2:],np.log(x[4:31]),1)[0] for x in D]
        i=np.where((md[:-1]>=pm)&(md[1:]<pm))[0][0];fr_=(np.log(md[i])-np.log(pm))/(np.log(md[i])-np.log(md[i+1]));el=float(np.exp(np.log(E[i])+fr_*(np.log(E[i+1])-np.log(E[i]))))
        Dx=Dl(sky(z,s,el,v));r[v]=dict(E_level=el,slope=float(np.polyfit(l[2:],np.log(Dx[4:31]),1)[0]),chi2=chi(Dx),r2=r2(Dx),D2_over_P=float(Dx[2]/Pl[0,1]),D30=float(Dx[30]),D4=float(Dx[4]))
    out[k]=r
json.dump(out,open(O+'flow_levelmatch.json','w'))
for v in ('grow','fixed'):
    for q in ('E_level','slope','chi2','r2','D2_over_P','D4','D30'):
        a=[o[v][q] for o in out.values()];print(v,q,np.percentile(a,[16,50,84]).round(4))
