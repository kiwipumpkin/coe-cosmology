# Phase4b (NEW): same model as phase4/k_scan.sky, but internal guards widened (Chris: no 2 Gpc limit).
# Old guards: 1<d<60, 1e-3<eta<8, 1e-4<mc<80, 1e-4<LD<3. New: 1e-3<=d<D/1.05, 1e-5<eta<1e3, 1e-6<mc<1e4, 1e-6<LD<100.
import sys,numpy as np;sys.path.insert(0,'/workspace/coepdf/joint_solve/phase4');import k_scan as K
cw=K.cw;TH=K.TH;AGE=K.AGE
def sky(D,d,eta,mc,k,LD=None):
    if not(D>d*1.05 and 1e-6<mc<1e4 and 1e-5<eta<1e3 and 1e-3<=d):return None
    if LD is None:
        try:L=cw.L_from_delta_eta(D,d,eta,0.031*k)
        except Exception:return None
        LD=L/D
    if not(1e-6<LD<100):return None
    try:rr=cw.pc.run('eta',dict(DR=D/cw.RV,LD=LD,eta=eta),mc,0.0)
    except Exception:return None
    af=rr['af']
    if not np.all(np.isfinite(af)) or rr['vrmin']<=0:return None
    return af,rr,LD
