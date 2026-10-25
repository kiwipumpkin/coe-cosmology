import numpy as np,healpy as hp,json,matplotlib;matplotlib.use('Agg');import matplotlib.pyplot as plt
T0=2.7255e6;P=np.loadtxt('planckTT.txt');l=P[:,0].astype(int);Cl=2*np.pi*P[:,1]/(l*(l+1))
pl={int(a):float(np.sqrt((2*a+1)*c/(4*np.pi))) for a,c in zip(l,Cl) if a<=12};pl[1]=3362/np.sqrt(3)
d=json.load(open('web_rec_D42.json'));rows=[]
print(f"{'l':>3} {'smooth(1e-5)':>12} {'web(1e-5)':>10} {'final(1e-5)':>11} {'final uK':>9} {'Planck uK':>9} {'Planck 1e-5':>11}")
for i in range(12):
    L=i+1;s,w,f=d['smooth']['per_l'][i],d['web_only']['per_l'][i],d['smooth_plus_web']['per_l'][i]
    rows.append(dict(l=L,smooth_1e5=s,web_1e5=w,final_1e5=f,final_uK=f*1e-5*T0,planck_uK=pl[L],planck_1e5=pl[L]/T0/1e-5))
    print(f"{L:3d} {s:12.3f} {w:10.3f} {f:11.3f} {f*27.255:9.1f} {pl[L]:9.1f} {pl[L]/27.255:11.3f}"+("  (both: before motion removal; final has it removed)" if L==1 else ""))
json.dump(rows,open('multipoles_vs_planck.json','w'),indent=1)
m=np.load('finalmap_D42.npy');w=np.asarray(hp.remove_dipole(np.load('webrec_D42.npy')))
fig=plt.figure(figsize=(16,13))
def mark():
    hp.graticule(30,alpha=.3)
    for (L,b),s,c in (((264.02,48.25),'Planck hot','r'),((84.02,-48.25),'Planck cold','b'),((180,-15),'CoE','m'),((0,15),'anti-CoE','g')):
        hp.projscatter(L,b,lonlat=True,marker='*',s=250,c=c,edgecolors='k');hp.projtext(L,b+6,s,lonlat=True,fontsize=12,weight='bold')
hp.mollview(m/1e-5,fig=fig.number,sub=(2,1,1),cmap='RdBu_r',unit='dT/T ×1e5',title=f'FINAL: smooth sequential sky (mass b, CatWISE, 1+z≈1100) − dipole (motion) + web at recombination contrast\nRMS {np.std(m)/1e-5:.0f}e-5 — dominated by the unremoved CoE-axis quadrupole');mark()
hp.mollview(w/1e-5,fig=fig.number,sub=(2,1,2),cmap='RdBu_r',min=-20,max=20,unit='dT/T ×1e5',title=f'Web contribution alone (seed-3, contrasts grown back to recombination): RMS {np.std(w)/1e-5:.1f}e-5 = {np.std(w)*T0:.0f} μK');mark()
plt.savefig('final_map.png',dpi=110,bbox_inches='tight')
