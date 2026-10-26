# Phase3 B: sourced H0 table (values verified on the web 2026-10-25), mean and band.
import json,numpy as np
S=[ # name, H0, sig_minus, sig_plus, ref, used_in_mean
("Planck 2018 TT,TE,EE+lowE+lensing (LCDM)",67.4,0.5,0.5,"Planck Collab. 2020, A&A 641, A6 (arXiv:1807.06209) abstract",1),
("ACT DR6 P-ACT-LB (CMB+lensing+DESI DR1)",68.22,0.36,0.36,"Louis et al. 2025, JCAP 11, 062 (arXiv:2503.14452)",0),
("DESI DR2 BAO + BBN",68.51,0.58,0.58,"DESI Collab. 2025, arXiv:2503.14738",1),
("CCHP TRGB HST+JWST (Freedman+25)",70.39,1.94,1.94,"Freedman et al. 2025, ApJ 985, 203; 1.22 stat,1.33 sys,0.70 SN in quadrature",1),
("SH0ES Cepheid+SN Ia (Riess+22)",73.04,1.04,1.04,"Riess et al. 2022, ApJL 934, L7",1),
("TDCOSMO 2025 lensing time delays (+Pantheon+ Om)",71.6,3.3,3.9,"TDCOSMO 2025, A&A (doi 10.1051/0004-6361/202555801, arXiv:2506.03023)",1),
("Megamaser Cosmology Project (Pesce+20)",73.9,3.0,3.0,"Pesce et al. 2020, ApJL 891, L1",1),
("GW170817 standard siren",70.0,8.0,12.0,"LIGO/Virgo et al. 2017, Nature 551, 85",1)]
H=np.array([s[1] for s in S]);sm=np.array([s[2] for s in S]);sp=np.array([s[3] for s in S]);use=np.array([s[5] for s in S],bool)
sg=0.5*(sm+sp);w=1/sg[use]**2
mean=float(np.sum(w*H[use])/w.sum());emean=float(w.sum()**-0.5)
prec=use&(sg<=2.0)
low=float(np.min(H[prec]-sm[prec]));high=float(np.max(H[prec]+sp[prec]))
out=dict(sources=[dict(name=s[0],H0=s[1],minus=s[2],plus=s[3],ref=s[4],in_mean=bool(s[5])) for s in S],
 H0_mean_invvar=mean,H0_mean_err=emean,H0_low=low,H0_high=high,midpoint=0.5*(low+high),
 rule="mean = inverse-variance (symmetrised sigma) over in_mean sources (ACT excluded: shares Planck data); band = [min(H-sig), max(H+sig)] over in_mean sources with sigma<=2 (Planck, DESI, CCHP, SH0ES)",
 chi2_scatter=float(np.sum(w*(H[use]-mean)**2)),n=int(use.sum()))
json.dump(out,open('h0_sources.json','w'),indent=1);print(json.dumps({k:v for k,v in out.items() if k!='sources'},indent=1))
