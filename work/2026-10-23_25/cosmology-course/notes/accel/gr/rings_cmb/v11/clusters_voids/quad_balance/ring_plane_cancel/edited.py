import sys,os,json,pickle,numpy as np
OUT='/workspace/clump_scan/sequential_open';sys.path.insert(0,OUT);sys.path.insert(0,'/workspace/clump_scan/sequential');H=os.getcwd();os.chdir(OUT)
import seq_core as S
os.chdir(H)
f=np.load(f'{OUT}/final_comp_s3.npz');X0=f['Xp'];r0=f['rp'];isV=f['kp']=='void';M0=f['Mp'];grp=f['grp'].astype(int)
NG=grp.max()+1;vb=np.zeros(NG,int);vb[grp[isV]]=np.where(isV)[0];c0=X0[vb];rhoc0=S.rho_gpc(c0);U0=X0-c0[grp]
st=pickle.load(open(f'{OUT}/ckpt/minedit_comp_comp_s3.pkl','rb'))
g=grp;s=st['A'][g]**(1/3);c=c0[g]+st['D'][g];X=c+s[:,None]*st['U'];M=M0*st['A'][g]*S.rho_gpc(c)/rhoc0[g];r=r0*s
T=S.R-f['w'];dW=S.pattern(X,M)-f['w'];res=T-dW
th=np.radians(np.arange(360)) if not hasattr(S,'TH') else S.TH
def mm(v):m=np.arange(1,61)[:,None];a=np.hypot(2*np.mean(v*np.cos(m*th),-1),2*np.mean(v*np.sin(m*th),-1));return a/np.sqrt(2)/S.FLOOR
moved=np.linalg.norm(X-X0,axis=1)>1e-9;resz=np.abs(M/M0-1)>1e-9
print('final total',S.rms(res)/S.FLOOR,'stored curve end',st['curve'][-1],'edited bodies',int((moved|resz).sum()),'units',len(st['edited']))
a=mm(res);print('final per-m',np.round(a[:12],3).tolist(),'m13-60',np.sqrt(np.sum(a[12:]**2)))
np.savez('final_s3_edited.npz',X=X,M=M,r=r,res=res,T=T)
json.dump(dict(final_total=float(S.rms(res)/S.FLOOR),per_m=a[:12].tolist(),m13_60=float(np.sqrt(np.sum(a[12:]**2))),edited_bodies=int((moved|resz).sum()),units=len(st['edited'])),open('paper_final.json','w'),indent=1)
