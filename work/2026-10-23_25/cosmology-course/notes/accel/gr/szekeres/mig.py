import sys,json,pickle,numpy as np; sys.path.insert(0,'/workspace/cosmology-course/notes/paper_repro/grfit')
import szk, toy
from core import ltb as LTBfit, shifted, chi, r2, TH1, RV
V={r['i']:r for r in json.load(open('../valley/valley_cmb.json'))}; P=json.load(open('../valley/valley_points.json'))
out=pickle.load(open("mig.pkl","rb"))
for i in (0,2):
    p=P[i]; fac,Dt=toy.make(p); L=fac.func(*fac.args[:3],V[i]['rmin']); a=1-p['kin']/1.0; rin=(p['D']-p['kin']*RV)/RV; rout=(p['D']+p['kout']*RV)/RV
    m0=shifted(LTBfit(p['D'],p['mcc'],p['kin'],p['kout']),p['l0']); print(i,'LTB fit code chi2 %.4f R2 %.4f'%(chi(m0),r2(m0)),'shell (ring units)',rin,rout,'D',Dt,flush=True)
    for C in (0.9,-0.9):
        S=szk.SZ(L,C,rin,min(rout,L.rg[-1]))
        m=shifted(szk.migkas_af(S,Dt,TH1),p['l0']); out[(i,C)]=dict(chi2=chi(m),R2=r2(m),sig0=S.sig0,ratio=S.ratio,af=m)
        print(i,C,'sig0 %.4f limit ratio %.3f chi2 %.4f R2 %.4f'%(S.sig0,S.ratio,chi(m),r2(m)),flush=True)
        pickle.dump(out,open('mig.pkl','wb'))
