# Cosmic rheology: modes of motion on top of the primary expansion

## What this section does

The special relativity (SR) plus central-mass treatment describes one motion: radial expansion of matter away from the centre of expansion (CoE), slowed by the central mass. Real matter is not a smooth shell. It is arranged into voids, filaments (threads) and clumps, and it can also shear and rotate. Here we treat the matter as a continuous medium, the way rheology treats a deforming material, and ask which secondary motions can sit on top of the primary expansion and what sky pattern each one would leave. The multipole residuals of the SR fit (Section X) are the fingerprint these motions would have to reproduce. No new fits are reported in this section.

## Decomposing the motion

At any point, the velocity field $\mathbf{v}(\mathbf{x})$ can be split into its primary expansion and a perturbation,

$$
\mathbf{v}=\mathbf{v}_0(r)\,\hat{\mathbf{r}}+\delta\mathbf{v}(\mathbf{x}),
$$

where $\mathbf{v}_0$ is the radial CoE flow used in the SR fit. The local velocity gradient separates uniquely into three parts:

$$
\frac{\partial v_i}{\partial x_j}=\frac{1}{3}\,\Theta\,\delta_{ij}+\sigma_{ij}+\omega_{ij}.
$$

Here $\Theta=\nabla\cdot\mathbf{v}$ is the expansion rate, $\sigma_{ij}$ is the trace-free symmetric shear, which stretches matter into filaments and sheets, and $\omega_{ij}$ is the antisymmetric vorticity, or local rotation. The primary CoE flow carries only $\Theta$. Voids, threads and clumps are made by $\sigma_{ij}$ together with spatial changes in $\Theta$, and rotation is carried by $\omega_{ij}$.

## The medium is dust, not an elastic solid

If the matter shell behaved like a constant-density elastic solid at the mean cosmic density, $\rho\approx8.5\times10^{-27}\ \mathrm{kg\,m^{-3}}$, the residual harmonics would be standing normal modes, like Earth's free oscillations. Their frequencies would be

$$
f_\ell\simeq\frac{\sqrt{\ell(\ell+1)}}{2\pi R}\sqrt{\frac{K}{\rho}},
$$

where $K$ is the bulk modulus. Even with a radiation-fluid sound speed $c/\sqrt3$, the $m=3$ mode on the $R=14.3$ Gpc scale would have a period of roughly 50 to 100 billion years, so the pattern would look frozen today. After recombination, though, the matter is pressureless dust with essentially no stiffness, so $K\to0$ and there are no restoring oscillations. The appropriate picture is therefore a velocity ripple. Matter flows slowly under its own gravity and that of the central mass, and the energy is kinetic plus gravitational potential energy from the bunching. Such gravitational flows grow monotonically rather than oscillating.

## Voids, threads and clumps

In a pressureless medium, matter drains out of underdense regions (voids) and collects along filaments and in clumps. Each produces a peculiar velocity $\delta\mathbf{v}$ on top of the CoE flow. Matter streams outward from void centres and inward toward filaments and clumps. Through the Doppler factor $D_{\rm motion}$, those velocities change the observed temperature by

$$
\frac{\Delta T}{T}\approx\frac{\delta\mathbf{v}\cdot\hat{\mathbf{n}}}{c}.
$$

A smaller contribution comes from the gravitational redshift of the potential wells and hills.

The dominant residual modes, $m=3$ to 5, have angular wavelengths $\lambda_m=2\pi R/m\approx18$–$30$ Gpc at $R=14.3$ Gpc, which is horizon scale. Treated as velocity ripples in a shell holding the baryon fraction $\Omega_s=0.049$ of the critical mass, they need about $7\times10^{64}$ J, roughly 0.002% of the shell's kinetic energy. In other words, the fingerprint asks for only a small departure from smooth expansion. The four-clump tidal sketch (Figs. 20 and 21) is one concrete example. Four equal masses of $7.7\times10^{51}$ kg each, at 24.5 Gpc and placed symmetrically about the shell axis, reproduce 88% of the post-SR residual variance. The $m=6$ and $m=1$ pattern left after them shows that the true distribution is less symmetric than four equal clumps.

## Differential rotation

If the shell rotates, and especially if its angular speed varies with radius or latitude, the vorticity term $\omega_{ij}$ is non-zero and the flow is no longer symmetric about the axis. In the Fourier analysis about $\phi=15.16^\circ$,

$$
R(\ell)=a_0+\sum_m\bigl[c_m\cos m(\ell-\phi)+s_m\sin m(\ell-\phi)\bigr],
$$

any rotation shows up as non-zero sine coefficients $s_m$. In the present SR fit, $|s_m|\sim10^{-8}$, so the residual is axisymmetric. The current model contains no rotation, and the residual itself does not require any. That makes rotation a clean, separate test. Adding the one-cycle phase factor of Migkas et al. raised the geometry fit from $R^2\approx0.55$ to $0.82$. A rotating model refit to the Migkas data and run through the same SR pipeline would predict specific $s_m$ that the residual could confirm or rule out.

## How each motion maps to the multipoles

| Motion on top of the primary expansion | Main sky signature |
|---|---|
| Bulk flow of the observer's region | $m=1$ dipole |
| Large-scale shear or tidal stretching | $m=2$ quadrupole |
| Horizon-scale voids, threads and clumps | $m=3$–6, cosine about the axis if placed symmetrically |
| Smaller structures | $m\ge7$, partly mixed with the point-emitter Doppler tail |
| Rigid or differential rotation | sine terms $s_m$ that break the axisymmetry |

## Caveats

This is a qualitative map, not a solved model. The kinetic-energy estimate assumes velocity ripples in a uniform shell with $\Omega_s=0.049$. The clump sketch assumes a steady tidal pull over 13.8 billion years, so its masses are uncertain by a factor of a few. Rotation has not been fitted. A full treatment would solve the continuity, Euler and Poisson equations for pressureless matter on the CoE background and compare the resulting $c_m$ and $s_m$ with the residual fingerprint.

## Suggested figures

- `notes/clump-fit-residual.png`: post-SR residual, four-clump model, and the leftover after the clumps.
- `notes/cmb-modes-sheet.png`: polar shapes of $m=1$–12, which show what each multipole looks like on the sky.
