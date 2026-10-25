#!/usr/bin/env python3
"""Reproducible wide/narrow ell=220 inverse-designed exterior-ring experiment.
Requires numpy, scipy, matplotlib.
Outputs two PNG figures and one comparison PNG in the current directory.

IMPORTANT: This computes gravitational-POTENTIAL angular power, NOT CMB temperature
anisotropy. It prescribes the target power and solves for signed shell contrasts.
Negative contrast means underdensity relative to a positive background.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import eval_legendre
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

L_MAX=420
N_RINGS=480
R=1.0
SOURCE_RADIUS=1.012*R
PEAK_L=220
WIDTHS=(32,8)
ell=np.arange(1,L_MAX+1)
mu,quad_w=np.polynomial.legendre.leggauss(N_RINGS)
P=np.array([eval_legendre(int(l),mu) for l in ell])

def solve(sigma):
    target=np.exp(-0.5*((ell-PEAK_L)/sigma)**2)
    target/=target.max()
    # Axisymmetric source ring expansion:
    # f(mu)=sum_l (2l+1)/2 * q_l * P_l(mu)
    # q_l=int_-1^1 f(mu)P_l(mu)dmu
    # Potential on R=1 due to exterior shell d>1:
    # Phi_l ∝ q_l/d^(l+1)
    # D_l ∝ l(l+1)(2l+1) |Phi_l|^2, up to common constants.
    q=np.sqrt(target/(ell*(ell+1)*(2*ell+1)))*SOURCE_RADIUS**(ell+1)
    contrast=((2*ell+1)/2*q)@P
    # Normalize density contrast only; normalization cancels in power shape.
    contrast/=np.max(np.abs(contrast))
    q_recovered=P@(quad_w*contrast)
    phi=q_recovered/SOURCE_RADIUS**(ell+1)
    power=ell*(ell+1)*(2*ell+1)*np.abs(phi)**2
    power/=power.max()
    return target,contrast,power

def plot_case(sigma, filename):
    target,contrast,power=solve(sigma)
    fig=plt.figure(figsize=(15,7.8),facecolor="white",layout="constrained")
    gs=fig.add_gridspec(1,3,width_ratios=[1.25,1.1,1.25])
    ax=fig.add_subplot(gs[0,0],projection="3d")
    az=np.linspace(0,2*np.pi,180)
    # All 480 physical latitude rings; slight line transparency to resolve overlaps.
    for j in range(N_RINGS):
        z=SOURCE_RADIUS*mu[j]
        radius=SOURCE_RADIUS*np.sqrt(1-mu[j]**2)
        color=plt.get_cmap("RdBu_r")((contrast[j]+1)/2)
        ax.plot(radius*np.cos(az),radius*np.sin(az),
                np.full_like(az,z),color=color,lw=.65,alpha=.75)
    ax.set(xlim=(-1.1,1.1),ylim=(-1.1,1.1),zlim=(-1.1,1.1))
    ax.set_box_aspect((1,1,1));ax.view_init(elev=24,azim=35)
    ax.set_title(f"480 latitude rings\nExterior shell r={SOURCE_RADIUS:.3f} R")
    ax.set_axis_off()

    ax2=fig.add_subplot(gs[0,1])
    # Meridian slice: each line is the cross section of an entire ring.
    for j in range(N_RINGS):
        z=SOURCE_RADIUS*mu[j]
        x=SOURCE_RADIUS*np.sqrt(1-mu[j]**2)
        color=plt.get_cmap("RdBu_r")((contrast[j]+1)/2)
        ax2.plot([-x,x],[z,z],color=color,lw=.8)
    ax2.add_patch(plt.Circle((0,0),R,fill=False,color="black",lw=1))
    ax2.set_aspect("equal");ax2.set(xlim=(-1.1,1.1),ylim=(-1.1,1.1))
    ax2.set_xlabel("x / R");ax2.set_ylabel("z / R")
    ax2.set_title("Meridional ring cross-sections")

    ax3=fig.add_subplot(gs[0,2])
    ax3.plot(ell,target,color="black",lw=2.2,label="Prescribed target")
    ax3.plot(ell,power,color="#d24a3b",ls="--",lw=1.8,label="Recovered from rings")
    ax3.axvline(PEAK_L,color=".5",ls=":")
    ax3.set(xlim=(100,340),ylim=(0,1.08),xlabel="Spherical multipole ℓ",
            ylabel="Normalized gravitational-potential Dℓ",
            title=f"Peak ℓ={PEAK_L}; σℓ={sigma}\nFWHM ≈ {2.355*sigma:.0f}")
    ax3.grid(alpha=.2);ax3.legend(fontsize=9)
    fig.suptitle("INVERSE-DESIGNED GRAVITATIONAL POTENTIAL — NOT CMB TEMPERATURE",fontsize=14)
    sm=ScalarMappable(norm=Normalize(-1,1),cmap="RdBu_r")
    cb=fig.colorbar(sm,ax=[ax,ax2],orientation="horizontal",shrink=.65,pad=.06)
    cb.set_label("Signed density contrast (normalized): blue = void, red = overdensity")
    fig.savefig(filename,dpi=180,bbox_inches="tight")
    plt.close(fig)
    return target,power

def main():
    outputs=[]
    fig,ax=plt.subplots(figsize=(10,5.5),layout="constrained")
    for sigma,color in zip(WIDTHS,("#6d6d6d","#c9473c")):
        path=f"rings_ell220_sigma{sigma}.png"
        target,power=plot_case(sigma,path)
        outputs.append(path)
        ax.plot(ell,power,color=color,lw=2,label=f"Recovered σℓ={sigma}")
    ax.axvline(PEAK_L,color="black",ls=":",lw=1)
    ax.set(xlim=(100,340),ylim=(0,1.08),xlabel="Spherical multipole ℓ",
           ylabel="Normalized gravitational-potential Dℓ",
           title="Wide versus narrow ℓ=220 potential-power peaks")
    ax.grid(alpha=.2);ax.legend()
    comparison="rings_ell220_spectrum_comparison.png"
    fig.savefig(comparison,dpi=180,bbox_inches="tight")
    plt.close(fig)
    outputs.append(comparison)
    print("Created:",*outputs,sep="\n  ")

if __name__=="__main__":
    main()
