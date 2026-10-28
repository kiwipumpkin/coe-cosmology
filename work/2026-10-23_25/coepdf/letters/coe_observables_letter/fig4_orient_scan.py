# Orientation scan. Convention = orient_lee/lee.py: R = rot(coe,roll) @ Rz(lc-180) @ rot(y,bc); model CoE axis (180,0) -> (lc,bc).
# Critical pts: model-frame b=0, l=180+-theta, theta=arccos(R/2D), plus antipodes. Offset = min angle to AoE.
# B planarity: B is axisymmetric about CoE (min-eig ratio <1e-11), so t_l about AoE depends only on angle(CoE,AoE); evaluated by measuring t_l of B about axes at that angle from its own axis. Roll does not enter.
import sys,json;sys.path.insert(0,'.');from common import *
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
exec(open('/workspace/coepdf/joint_solve/phase5/orient_lee/lee.py').read().split("if __name__")[0].split("MM=json")[0])
def vec(l,b):l,b=np.radians(l),np.radians(b);return np.stack([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)],-1)
def rot(ax,a):
  ax=ax/np.linalg.norm(ax);a=np.radians(a);K=np.array([[0,-ax[2],ax[1]],[ax[2],0,-ax[0]],[-ax[1],ax[0],0]]);return np.eye(3)+np.sin(a)*K+(1-np.cos(a))*K@K
A=vec(260.,60.);Z=np.array([0,0,1.]);Y=np.array([0,1.,0])
def off(lc,bc,roll,DR):
  R1=rot(Z,lc-180)@rot(Y,bc);c=R1@vec(180,0);Rm=rot(c,roll)@R1
  th=np.degrees(np.arccos(1/(2*DR)));P=vec(np.array([180+th,180-th]),np.zeros(2))@Rm.T
  return float(np.degrees(np.arccos(min(1,np.abs(P@A).max()))))
DRS={'phase4b D/R=27.34/21.74':27.34/21.74,'figure D/R=5':5.0}
rolls=np.arange(0,360.1,1.);LC=np.arange(150,210.1,2);BC=np.arange(-45,45.1,2)
assert abs(np.degrees(np.arccos(rot(Y,-15)@vec(180,0)@vec(180,-15)))) < 1e-6
# B planarity vs angle to AoE
alB=ALM['CoE_bg_only'];perp=np.cross(nco,Z);perp/=np.linalg.norm(perp);angs=np.arange(0,180.1,2)
g={l:[] for l in (2,3)}
for a in angs:
  n=np.cos(np.radians(a))*nco+np.sin(np.radians(a))*perp;ax=lb(n) if True else None
  lo,la=hp.vec2ang(n,lonlat=True)
  for l in (2,3):g[l].append(tl_about(alB,(float(lo[0]),float(la[0])),l)[0])
def tB(lc,bc,l):a=np.degrees(np.arccos(np.clip(vec(lc,bc)@A,-1,1)));return np.interp(a,angs,g[l])
out={'convention':'lee.py; paper orientation = (180,-15), roll 45','DR':DRS}
fig,ax=plt.subplots(3,3,figsize=(18,15));ax=ax.ravel()
def band(a):a.axvspan(45,80,color='y',alpha=.3,label='expected roll 45-80');a.axvline(60,color='r',ls='--',label='roll 60')
for i,(nm,DR) in enumerate(DRS.items()):
  a=ax[i]
  for (lc,bc) in [(180,-15),(180,0),(180,15)]:
    o=[off(lc,bc,r,DR) for r in rolls];a.plot(rolls,o,label=f'CoE ({lc},{bc})');out[f'{nm} offset vs roll CoE({lc},{bc})']={'roll60':o[60],'roll300':o[300],'min_all':[min(o),float(rolls[int(np.argmin(o))])],'min_45_80':[min(o[45:81]),float(rolls[45+int(np.argmin(o[45:81]))])]}
  band(a);a.set_xlabel('roll (deg)');a.set_ylabel('AoE offset of nearest critical pt (deg)');a.set_title(nm);a.set_xlim(0,360);a.legend(fontsize=7)
  for j,(ttl,fn) in enumerate([('roll = 60',lambda lc,bc:off(lc,bc,60,DR)),('min over roll 45-80',lambda lc,bc:min(off(lc,bc,r,DR) for r in range(45,81))),('min over roll 0-360',lambda lc,bc:min(off(lc,bc,r,DR) for r in range(0,360,2)))]):
    if i==1 and j>=1:continue  # slot ax[7] is redrawn as the t_3 map below; skip so no stray colorbar
    G=np.array([[fn(lc,bc) for lc in LC] for bc in BC]);a=ax[3+3*i+j] if j<2 else ax[2]
    im=a.pcolormesh(LC,BC,G,shading='nearest',cmap='viridis_r');plt.colorbar(im,ax=a,label='AoE offset (deg)')
    a.contour(LC,BC,G,[1,5,10],colors='w',linewidths=.8);a.plot(180,-15,'r*',ms=12);a.plot(180,0,'w+',ms=12)
    a.set_xlabel('CoE l (deg)');a.set_ylabel('CoE b (deg)');a.set_title(f'{nm}: {ttl}');a.invert_xaxis()
    k=np.unravel_index(np.argmin(G),G.shape);out[f'{nm} map {ttl}']={'min':float(G[k]),'at_lc_bc':[float(LC[k[1]]),float(BC[k[0]])],'at(180,-15)':float(G[list(BC).index(-15),list(LC).index(180)]),'at(180,-15)':float(G[list(BC).index(-15),list(LC).index(180)])}
# B planarity panels
a=ax[5];a.cla()
for l in (2,3):a.plot(rolls,[tB(180,-15,l)]*len(rolls),label=f't_{l} CoE bg, CoE (180,-15)')
band(a);a.set_xlim(0,360);a.set_ylim(0,1);a.set_xlabel('roll (deg)');a.set_ylabel('t_l about AoE');a.set_title('CoE bg planarity vs roll (axisymmetric: roll-independent)');a.legend(fontsize=7)
for l,a in [(2,ax[8]),(3,ax[7])]:
  a.cla();G=np.array([[tB(lc,bc,l) for lc in LC] for bc in BC]);im=a.pcolormesh(LC,BC,G,shading='nearest',cmap='magma');plt.colorbar(im,ax=a,label=f't_{l} about AoE')
  a.plot(180,-15,'c*',ms=12);a.set_xlabel('CoE l');a.set_ylabel('CoE b');a.set_title(f'CoE bg t_{l} about AoE vs CoE direction (any roll)');a.invert_xaxis()
  out[f'B t_{l} vs CoE dir']={'at(180,-15)':float(tB(180,-15,l)),'max_in_grid':float(G.max()),'min_in_grid':float(G.min())}
out['B t_l vs angle(CoE,AoE)']={str(l):dict(zip(map(str,angs[::5]),np.round(g[l][::5],3).tolist())) for l in (2,3)}
plt.tight_layout();plt.savefig('fig4_orientation_scan.png',dpi=100);json.dump(out,open('fig4_values.json','w'),indent=1);print(json.dumps(out,indent=1))
