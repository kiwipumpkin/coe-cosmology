"""Redshift decomposition along past light cone: ln(1+z) = -ln[g_o(1-b_o.k)] + ln[g_e(1-b_e.k)] + (Phi_o-Phi_e)/c^2 (Newtonian shared)."""
import os,sys,json,numpy as np,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
law,mode=sys.argv[1],sys.argv[2]; os.environ['VMODE']='gtc'; sys.argv=['x',law,mode]
exec(open('fit_lc.py').read().split("ck=f'fwd")[0])
x=np.array(json.load(open(f'fwd/lc_{law}_{mode}_gtc_ckpt.json'))['x']); P,o,af,chi,l0=evl(x)
io=o['io']; D=o['r'][io]; Uo=o['U'][io]; tn=o['t']; uo=Uo/F.C; bo=uo/np.sqrt(1+uo**2); go=np.sqrt(1+uo**2); pho=F._phi(o,o['r'])[io]
ZT=[0.01,0.1,0.3,1,10,100,1100]; Dl=np.geomspace(1e-6,0.999,600)*tn; out={}
for thd in (0,90,180):
    th=np.radians(thd); nx,ny=np.cos(th),np.sin(th); kx,ky=-nx,-ny; rows=[]
    for dl in Dl:
        r,U=F._state_at(o,tn-dl); live=r>1.1e-6*F.R0; x_=D+F.C*dl*nx; y_=F.C*dl*ny; re=np.hypot(x_,y_)
        o2=np.argsort(r[live]); rr=r[live][o2]
        if re>rr.max() or re<rr.min(): rows.append([np.nan]*6); continue
        ue=np.interp(re,rr,U[live][o2])/F.C; ph=np.interp(re,rr,F._phi(o,r)[live][o2]); Me=np.interp(re,rr,o['GM'][live][o2])
        be=ue/np.sqrt(1+ue**2); ge=np.sqrt(1+ue**2)
        To=-np.log(go*(1-bo*kx)); Te=np.log(ge*(1-be*(x_*kx+y_*ky)/re)); Tg=(pho-ph)/F.C**2
        rows.append([To,Te,Tg,To+Te+Tg,re/D,Me/o['GM'][io]])
    R=np.array(rows); lz=R[:,3]; res=[]
    for z in ZT:
        k=np.where(np.isfinite(lz)&(lz>=np.log1p(z)))[0]
        if not len(k): res.append(dict(z=z,reached=False)); continue
        k=k[0]; a=R[k]; res.append(dict(z=z,reached=True,lookback_frac=float(Dl[k]/tn),obs=float(a[0]),emit=float(a[1]),grav=float(a[2]),
            f_obs=float(a[0]/a[3]),f_emit=float(a[1]/a[3]),f_grav=float(a[2]/a[3]),r_emit_over_D=float(a[4]),Menc_emit_over_obs=float(a[5]),
            min_r_on_path_over_D=float(np.nanmin(R[:k+1,4]))))
    out[thd]=res
# anisotropy decomposition (0 vs 180 deg difference of each term) at each z
an=[]
for i,z in enumerate(ZT):
    a,b=out[0][i],out[180][i]
    if a.get('reached') and b.get('reached'): an.append(dict(z=z,d_obs=a['obs']-b['obs'],d_emit=a['emit']-b['emit'],d_grav=a['grav']-b['grav']))
json.dump(dict(law=law,mode=mode,U_o_c=uo,dirs=out,aniso=an),open(f'fwd/zdecomp_{law}_{mode}.json','w'),indent=1,default=float)
fig,ax=plt.subplots(1,2,figsize=(16,8))
for thd,m in ((0,'o-'),(180,'s--')):
    rr=[r for r in out[thd] if r.get('reached')]; z=[r['z'] for r in rr]
    for k,c in (('f_obs','C0'),('f_emit','C1'),('f_grav','C2')): ax[0].semilogx(z,[r[k] for r in rr],m,c=c,lw=3,ms=9,label=f'{k[2:]} θ={thd}')
    ax[1].loglog(z,[r['r_emit_over_D'] for r in rr],m,lw=3,ms=9,label=f'r_emit/D θ={thd}')
ax[0].set(xlabel='z',ylabel='fraction of ln(1+z)');ax[1].set(xlabel='z',ylabel='emitter distance from CoE / D')
for a in ax: a.legend(fontsize=13);a.grid();a.xaxis.label.set_size(20);a.yaxis.label.set_size(20);a.tick_params(labelsize=15)
fig.suptitle(f'Redshift decomposition, {law} {mode}, U_o/c={uo:.3g}',fontsize=20);fig.tight_layout();fig.savefig(f'fwd/zdecomp_{law}_{mode}.png',dpi=80)
print(json.dumps(an,default=float)[:800])
