import os,json,numpy as np,warnings;warnings.filterwarnings('ignore')
import matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
H=os.getcwd();os.chdir('/workspace/cosmology-course/notes/accel/gr/three_vel/density')
exec(open('/workspace/cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/migkas_prof.py').read().split('rng=np.random')[0]);os.chdir(H)
R=[d for d in json.load(open('fitDM_results.json')) if 'chi2' in d];P=np.loadtxt('../../../reproduce_cmb/planckTT.txt')
DS=sorted({d['DRV'] for d in R});FS=sorted({d['frac'] for d in R})
fig,ax=plt.subplots(2,3,figsize=(16,9))
for i,sh in enumerate(['red','ref']):
    for j,(key,lab) in enumerate([('chi2','Migkas chi2 (l0 free)'),('D2','log10 D2 [uK^2]'),('pschi2','log10 PS chi2 l=2..30')]):
        M=np.full((len(DS),len(FS)),np.nan)
        for d in R:
            if d['shape']==sh:
                v=d['chi2'] if key=='chi2' else np.log10(d['Dl'][2] if key=='D2' else d['pschi2']);M[DS.index(d['DRV']),FS.index(d['frac'])]=v
        a=ax[i,j];im=a.imshow(M,origin='lower',aspect='auto',cmap='viridis');plt.colorbar(im,ax=a)
        for y in range(len(DS)):
            for x in range(len(FS)):
                if np.isfinite(M[y,x]):a.text(x,y,'%.2f'%M[y,x],ha='center',va='center',fontsize=6,color='w')
        a.set_xticks(range(len(FS)));a.set_xticklabels(['%g'%f for f in FS],rotation=45,fontsize=7);a.set_yticks(range(len(DS)));a.set_yticklabels(['%g'%v for v in DS])
        a.set_xlabel('frac');a.set_ylabel('D/R (DRV)');a.set_title('%s: %s (blank=failed)'%(sh.upper(),lab),fontsize=9)
plt.tight_layout();plt.savefig('fitDM_maps.png',dpi=80);plt.close()
l=np.arange(31)
for sh in ['red','ref']:
    SS=sorted([d for d in R if d['shape']==sh],key=lambda d:(d['DRV'],d['frac']));plt.figure(figsize=(9,6));cm=plt.cm.viridis
    for k,d in enumerate(SS):plt.plot(l[2:],np.array(d['Dl'])[2:],color=cm(k/max(1,len(SS)-1)),lw=0.8,label='D/R %g f %g'%(d['DRV'],d['frac']))
    plt.errorbar(P[:29,0],P[:29,1],yerr=[P[:29,2],P[:29,3]],fmt='ko',ms=3,label='Planck 2018')
    plt.yscale('log');plt.xlabel('l');plt.ylabel('D_l [uK^2]');plt.title('%s shape: D_l, all %d scan runs'%(sh.upper(),len(SS)));plt.legend(fontsize=4,ncol=4);plt.savefig('fitDM_Dl_overlay_%s.png'%sh,dpi=80);plt.close()
def best(sh,key):return min([d for d in R if d['shape']==sh],key=lambda d:d[key])
Ls=np.arange(0,361,1.0);Lsave=L.copy()
for sh in ['red','ref']:
    for key in ['chi2','pschi2']:
        d=best(sh,key);tag='%s_best_%s'%(sh,'migkas' if key=='chi2' else 'pschi2');fn='fitDM_runs/%s_D%g_f%g.npz'%(sh,d['DRV'],d['frac'])
        fig,ax=plt.subplots(3,1,figsize=(8,12));th=np.array(d['theta'])
        L=np.mod(Ls,360);L[-1]=360.;m=model(th,d['ratio']);L=Lsave
        ax[0].errorbar(L,Y,yerr=S,fmt='ko',ms=3,label='71 Migkas bins');ax[0].plot(Ls,m,'r',label='model chi2 %.2f R2 %.3f l0 %.1f'%(d['chi2'],d['R2'],d['l0']));ax[0].set_xlabel('gal. l [deg]');ax[0].set_ylabel('A/A_all');ax[0].legend(fontsize=7)
        o=np.load(fn)['o'];o=o[np.argsort(o[:,0])];cs=CubicSpline(o[:,0],np.log(o[:,1]));ps=np.linspace(0,180,721);y=np.exp(-cs(ps));mu=np.cos(np.radians(ps))
        w=np.sin(np.radians(ps));y0=y/np.trapezoid(y*w,ps)*np.trapezoid(w,ps)-1;a1=1.5*np.trapezoid(y0*mu*w,ps)/np.trapezoid(w,ps)*2;y1=y0-a1*mu
        ax[1].plot(ps,y0,label='monopole removed');ax[1].plot(ps,y1,label='monopole+dipole removed');ax[1].set_xlabel('psi from CoE [deg] (CoE at 0)');ax[1].set_ylabel('dT/T');ax[1].legend(fontsize=7)
        ax[2].plot(l[2:],np.array(d['Dl'])[2:],'r-o',ms=3,label='model, PS chi2 %.3g'%d['pschi2']);ax[2].errorbar(P[:29,0],P[:29,1],yerr=[P[:29,2],P[:29,3]],fmt='ko',ms=3,label='Planck 2018');ax[2].set_yscale('log');ax[2].set_xlabel('l');ax[2].set_ylabel('D_l [uK^2]');ax[2].legend(fontsize=7)
        fig.suptitle('%s  D/R=%g frac=%g ratio=%.6f'%(tag,d['DRV'],d['frac'],d['ratio']));plt.tight_layout();plt.savefig('fitDM_%s.png'%tag,dpi=80);plt.close()
        print(tag,d['DRV'],d['frac'],d['chi2'],d['pschi2'])
