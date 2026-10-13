import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
rng=np.random.default_rng(3); plt.rcParams.update({'font.size':20})
HA=np.radians(25); R0,R1=0.12,1.0; N=600
xs=np.linspace(-0.45,0.45,N); ys=np.linspace(0,1.02,N); h=xs[1]-xs[0]
X,Y=np.meshgrid(xs,ys); RR=np.hypot(X,Y); TT=np.arctan2(X,Y)
mask=(RR>R0)&(RR<R1)&(abs(TT)<HA)
# Lagrangian particles, polar, uniform area density
M=2500000; qr=np.sqrt(rng.uniform(R0**2,R1**2,M)); qt=rng.uniform(-HA,HA,M); w=qr*0+1
Nsh=14; dsh=(R1-R0)/Nsh; shell=np.floor((qr-R0)/dsh); phase_s=rng.uniform(0,2*np.pi,Nsh+1)
lam=3*dsh*0.35; k=2*np.pi/lam  # wavelength ~3x shell thickness (thickness 0.35*spacing)
s=qr*qt
growth=np.clip((qr-R0)/(R1-R0),0,1)**1.2  # inner shells stay homogeneous
def dep(rr,tt,wt,sm=1.2):
    x=rr*np.sin(tt); y=rr*np.cos(tt)
    H,_,_=np.histogram2d(y,x,bins=[N,N],range=[[ys[0]-h/2,ys[-1]+h/2],[xs[0]-h/2,xs[-1]+h/2]],weights=wt)
    H=gaussian_filter(H,sm); H/=H[mask].mean(); return H
def frame(rr,tt,title,tl,fn,sm=1.2,vmax=3):
    D=dep(rr,tt,w,sm); D=np.where(mask,D,np.nan)
    f,a=plt.subplots(figsize=(9,10)); im=a.pcolormesh(X,Y,D,cmap='inferno',vmin=0,vmax=vmax,shading='auto')
    a.set_aspect('equal'); a.axis('off'); a.plot(0,0,'c^',ms=14); a.text(0.02,-0.01,'CoE',color='c',fontsize=18,va='top')
    a.set_title(title+'\n'+tl,fontsize=24); c=f.colorbar(im,ax=a,fraction=0.04,pad=0.02); c.set_label('density / mean',fontsize=20)
    a.text(0.5,0.0,'illustrative toy model, schematic units',ha='center',fontsize=14,color='gray',transform=a.transAxes)
    f.tight_layout(); f.savefig(fn,dpi=110); plt.close(f); return fn
# shell compression: particles bunch into thin shells
cen=R0+(shell+0.5)*dsh
def shells(th): return cen+(qr-cen)*th
A=0.9  # tangential compression amplitude -> density 0.1..1.9 x mean
fr=[]
fr.append(frame(shells(0.35),qt,'1. Flat onion shells','t ≈ recombination','f1.png'))
# ripples: radial wrinkle + tangential bunching
rip=np.sin(k*s+phase_s[shell.astype(int)])
tan=lambda a: qt+(a*growth*np.cos(k*s+phase_s[shell.astype(int)])/k)/qr
fr.append(frame(shells(0.35)+0.25*dsh*growth*rip,tan(A),'2. Shells ripple (Vishniac wrinkles)','t₁: λ ≈ 3× shell thickness','f2.png'))
# damping: radial ripple decays, shells thicken/merge, tangential pattern frozen
fr.append(frame(shells(0.75)+0.06*dsh*growth*rip,tan(A),'3. Ripples damp, pattern freezes','t₂','f3.png'))
fr.append(frame(qr,tan(A),'4. Lattice of cluster & void seeds','t₃: shells merged','f4.png',sm=2.5))
# Zel'dovich growth from seed field + small random modes
x0=qr*np.sin(tan(A)); y0=qr*np.cos(tan(A))
d0=(dep(qr,tan(A),w,4)-1)*0.5; d0=np.where(mask,d0,0)+1.2*gaussian_filter(rng.normal(size=(N,N)),10)/np.std(gaussian_filter(rng.normal(size=(N,N)),10))*0.15*np.where(mask,1,0)*(RR/R1)
K=np.fft.fftfreq(N,h)*2*np.pi; KX,KY=np.meshgrid(K,K); K2=KX**2+KY**2; K2[0,0]=1
dk=np.fft.fft2(d0); px=np.real(np.fft.ifft2(1j*KX*dk/K2)); py=np.real(np.fft.ifft2(1j*KY*dk/K2))
ix=np.clip(((x0-xs[0])/h).astype(int),0,N-1); iy=np.clip(((y0-ys[0])/h).astype(int),0,N-1)
psx,psy=px[iy,ix],py[iy,ix]; sc=0.010/np.std(np.hypot(psx,psy))
for i,(Dg,lab) in enumerate([(0.6,'5. Voids expand, matter drains to walls'),(1.1,'6. Walls collapse into filaments'),(1.7,'7. Cosmic web: nodes at intersections')]):
    xx=x0+Dg*sc*psx; yy=y0+Dg*sc*psy
    fr.append(frame(np.hypot(xx,yy),np.arctan2(xx,yy),lab,f't{4+i}: growth factor {Dg}',f'f{5+i}.png',sm=1.5,vmax=3))
import PIL.Image as I
ims=[I.open(f) for f in fr]; W,Hh=ims[0].size; g=I.new('RGB',(W*4,Hh*2),'white')
for i,im in enumerate(ims): g.paste(im,((i%4)*W,(i//4)*Hh))
g.save('filmstrip.png'); ims[0].save('evolution.gif',save_all=True,append_images=ims[1:],duration=1500,loop=0)
