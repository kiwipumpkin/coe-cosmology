import os,re,sys
sys.path.insert(0,'.');from cap import captions
from PIL import Image
W='/workspace/'
FIG='/workspace/coepdf/fig/'
capnow=captions('/workspace/coepdf/main_rev5.tex')
capmorn=captions('/workspace/coepdf/chk/main_rev5_pre_fullrewrite.tex')
MK=W+'figs_letter_exact_python/makers/'
SF=W+'redweb2/stochfield/'
DP=W+'coepdf/dens_plots.py'
CP='copy time onto the box (true original likely earlier)'
# order, name, files, date, source, morning, now, scripts, notes
F=[
('Passing-universe domains',['f05.png'],'Wed Oct 21 2026, 4:19 PM CT','git commit ae85632 (coe-cosmology-gh), adding mc_domain_overlay_lin.png',18,18,[W+'redweb2/realweb/mc/mc_overlay.py'],
 'Domains at R = 4, 8, 14 Gpc overlaid, with spreads over Monte-Carlo realisations of the realistic web. Derived from the redweb2/realweb Monte-Carlo runs; no earlier figure in this book is an input. Script saves mc_domain_overlay_lin.png, which is md5-identical to f05.png.'),
('First stochastic-field fit (stage A)',['f15.png'],'Thu Oct 22 2026, 12:42 AM CT','git commit c18eae2, adding stochfield/a_cl_l30_bg.png',12,21,[W+'coepdf/book/gitsrc/figs.py',W+'coepdf/book/gitsrc/fwd.py'],
 'Gaussian random field riding the CoE flow, R = 4 Gpc, fitted to Planck low multipoles. Lives in the stochfield/ work; builds on the CoE flow geometry. Script (from GitHub kiwipumpkin/coe-cosmology, commit c18eae2, work/2026-10-21_22/stochfield/figs.py) saves a_cl_{tag}.png; run with tag l30_bg it writes a_cl_l30_bg.png, which is md5-identical to f15.png. fwd.py is its imported helper.'),
('Parameter degeneracies',['f14.png'],'Thu Oct 22 2026, 12:42 AM CT','git commit c18eae2 (same commit as #2; ordered after it by paper number), adding stochfield/d_degeneracy_l30_bg.png',22,23,[SF+'deg.py'],
 'Delta-chi2 contours of shell position vs amplitude for the stochastic-field fit. Depends on the stage-A fit (#2). Script writes d_degeneracy_l30_bg.png (md5 match to f14.png).'),
('Same source at R = 4, 5, 6 Gpc',['f03.png'],'Thu Oct 22 2026, 12:42 AM CT','git commit c18eae2, adding stochfield/acoustic_ladder_456.png (pixel-identical to f03.png)',21,22,[W+'coepdf/book/gitsrc/acoustic_ladder456.py',W+'coepdf/book/gitsrc/acoustic_model.py'],'Peak positions vs shell radius for the same source parameters. Script (commit c18eae2, work/2026-10-21_22/stochfield/acoustic_ladder456.py) ends with savefig(\'acoustic_ladder_456.png\'); that committed PNG is pixel-for-pixel identical to f03.png (f03 is a re-saved copy). acoustic_model.py is its imported module.'),
('Low multipoles: LCDM, Planck, CoE one-sky',['f17.png'],'Thu Oct 22 2026, 12:42 AM CT','git commit c18eae2, adding stochfield/combo10.png (pixel-identical to f17.png)',20,27,[W+'coepdf/book/gitsrc/combo10.py'],'CoE one-sky realisation in its cosmic-variance band. Script (commit c18eae2, work/2026-10-21_22/stochfield/combo10.py) ends with savefig(\'combo10.png\'); that committed PNG is pixel-for-pixel identical to f17.png.'),
('Damped-harmonic fits to eps_l',['f18.png'],'Thu Oct 22 2026, 12:42 AM CT','git commit c18eae2, adding stochfield/resid_damped_7to9.png (pixel-identical to f18.png)',23,24,[],'Fits with 7, 8, 9 harmonics to eps_l. Original name resid_damped_7to9.png found in git; no .py in any commit or under /workspace contains resid_damped, so its script was never saved (likely run inline).'),
('Final CMB model vs Planck TT',['f01.png'],'Thu Oct 22 2026, 11:14 AM CT','git commit 2b5677e, adding acoustic_half_cl_fixed.png',24,25,[SF+'half_figs_fixed.py',SF+'half_fit.py',SF+'acoustic_fit.py',SF+'acoustic_model.py'],
 'CoE acoustic field plus damped half-scale term vs Planck TT, exact projection to l = 2500. Depends on the stochastic-field / acoustic fits (#2, #3). Helper modules half_fit, acoustic_fit and acoustic_model are included below.'),
('Density and t vs a (radiation + matter)',['d4_radmat.png'],'Fri Oct 23 2026, 12:37 AM CT','git commit 54f717f, adding math_pdf/fig/d4_radmat.png',3,3,[DP],
 'Eq. 1 density and Eq. 2 t(a) for a_eq = 1/3400, t0 = 13.8 Gyr. Analytic; no figure inputs. The four d*.png figures (#5–#8) come from one commit; ordered by paper number.'),
('CatWISE-locked density laws',['d1_migkas_laws.png'],'Fri Oct 23 2026, 12:37 AM CT','git commit 54f717f',5,5,[DP],
 'rho/rho_crit vs distance from the CoE with zoom on D +/- d, from the CatWISE-locked level-set laws. Input: CatWISE 3% dipole lock.'),
('Cosmic-web ingredients',['d2_web.png'],'Fri Oct 23 2026, 12:37 AM CT','git commit 54f717f',13,7,[DP],
 'Halo mass function, void size function and compensated void mass placement, recomputed from cosmo.py. Input to the realistic web (#14) and later web figures.'),
('Fitted P_delta terms and window W(chi)',['d3_stoch.png'],'Fri Oct 23 2026, 12:37 AM CT','git commit 54f717f',19,26,[DP],
 'Fitted P_delta split into its terms and the radial window. Depends on the stochastic-field fit (#2).'),
('Observables: Migkas A/A_all with bars and rotation overlay',['migkas_bars_rotation_overlay.png'],'Fri Oct 23 2026, 8:21 AM CT',CP,4,4,[W+'cosmology/figures/migkas_bars_rotation_overlay.py'],
 'Migkas A/A_all vs Galactic longitude with bars and rotation overlay. Input: Migkas cluster data (reads a spreadsheet via openpyxl).'),
('Sky-orientation sequence (four rows)',['Fig05_sky_orientation_four_row.jpg'],'Fri Oct 23 2026, 8:23 AM CT (8:23:32)',CP,7,10,[MK+'make_fig05_sky_orientation.py',MK+'paths.py'],
 'Untilted, +15 deg tilt, radial gradient, +80.8 deg roll. The maker stacks panels produced by four ring scripts it loads dynamically from /workspace/cosmology/figures/: make_coe_rings_notilt.py, make_coe_rings_migkas_pin_axis_matched.py, make_coe_rings_migkas_pin_axis_matched_grad.py, make_coe_rings_migkas_pin_axis_matched_grad_rot_fit.py (named, not reproduced). Caption says restored from the original preprint.'),
('Full Migkas A/A_all with model and residual bars (key plot)',['Fig04_amp_full.png'],'Fri Oct 23 2026, 8:23 AM CT (8:23:33)',CP+'; dated from an earlier same-name version',6,9,[MK+'make_fig04_amp_full.py',MK+'paths.py',MK+'physics.py'],
 'Full Migkas profile with model curves and residual bars. Input: Migkas data, CoE model from physics.py; follows the observables plot (#9).'),
('Geometry set-up',['Fig01a_geometry_top_test.png','Fig01a_geometry_bot_test.png','Fig01b_geometry.png'],'Fri Oct 23 2026, 8:23 AM CT (8:23:52)',CP+'; dated from an earlier same-name Fig01b_geometry.png',1,1,[W+'coepdf/scripts/fig01b_geometry.py'],
 'CoE, observer (Earth) and galaxy velocity vectors; D from CoE to observer. Top and middle panels: generating script not found (taken from the original preprint docx). Bottom panel (Fig01b) was redrawn on Oct 24-25 by the script below, so the caption phrase "Restored from the original preprint" no longer holds for it.'),
('Realistic cosmic web slab',['f06.png'],'Fri Oct 23 2026, 8:29 AM CT',CP,16,16,[],
 'Slab |z|<0.4 Gpc of a 20 Gpc box, seed 2000, 3904 clusters. Built from web ingredients (#7). Candidates searched: redweb2/stack_s2000.py writes red_web_s2000_stack.png (a spectrum stack, not this slab), so not used.'),
('Migkas 71 bins with CoE direct fits',['f07.jpg'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.584)',CP,15,12,[W+'cosmology-course/notes/accel/gr/rings_cmb/v11/clusters_voids/quad_balance/mass_family_scan/constrained/poly_plots.py'],
 'Red polynomial uses the same c_k as m(r) (#27). Caveat: poly_plots.py saves poly_plots.png / poly_spectrum.png, not f07.jpg directly; f07.jpg is taken to be its output per the inventory.'),
('Volume vs time (radiation + matter)',['f02.png'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.612)',CP+'; same second as #16, ordered by paper number',2,2,[],
 'Volumes of the D sphere, D+R sphere and outer sphere from a(t) of Eq. 2; follows from #5.'),
('Target profiles matched numerically',['f12.png'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.612)',CP,10,19,[],
 'Fixed mass at 1.55R plus fitted clumps/voids. Per Chris, matched to a different power spectrum from the web/ladder track.'),
('Axisymmetric harmonics on a ring',['f11.png'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.824)',CP,9,8,[],'Single harmonics l = 0..7 as density on a ring with the power bar at that l.'),
('Target peak matched numerically',['f13.png'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.880)',CP,11,20,[],'Non-axisymmetric source shell at r = 1.012R. Per Chris, made because the spectrum from the target-profiles fit (#16) was not there yet.'),
('Density modulation m(r)',['d5_mod.png'],'Fri Oct 23 2026, 8:30 AM CT (8:30:28.880)',CP+'; same instant as #21, ordered by paper number; redrawn since (last Oct 25 12:27 AM CT)',14,6,[W+'coepdf/d5_mod_earth.py',W+'coepdf/scripts/shp.py'],
 'm(r) with Earth at 42 Gpc, band 27.7-56.3 Gpc; illustrative unconstrained falloffs beyond 56.3 Gpc, cut at 80 Gpc. Original version came from coepdf/dens_plots.py (shown under #5). c_k shared with #14.'),
('A/B running-mean residuals and sky maps',['sky_angle.png','sky_mapA_al.png','sky_mapB_al.png'],'Sat Oct 24 2026, 4:40 PM CT','file creation, coepdf/fig/sky_angle.png',27,15,[],
 'Residuals vs model-frame angle and residual sky maps for A and B. Depends on the model fits (#11). coepdf/conv.py only references these files in LaTeX; it does not make them.'),
('Spot-scores table',['sky_tab_only.png'],'Sat Oct 24 2026, 4:46 PM CT','file creation, coepdf/fig/sky_tab_only.png',26,14,[],'Spot scores by sign of the tied residual (A vs B spots). Built from #23.'),
('Predicted MAX/MIN bands on Migkas Fig 8',['Fig06_overlay_one.png'],'Sat Oct 24 2026, 8:17 PM CT','file creation, coepdf/fig/Fig06_overlay_one.png',8,11,[],
 'Bands overlaid on Migkas Fig 8. Related makers exist (figs_letter_exact_python/makers/make_fig06_overlay_stack.py, make_fig06_matplotlib_clean.py) but write Fig06_overlay_stack.png / Fig06_matplotlib_clean.png, not this file.'),
('Web smoothed ladder (R = 0.5-14 Gpc)',['f19mc3b_stack.png'],'Sat Oct 24 2026, 11:29 PM CT','file creation, coepdf/fig/f19mc3b_stack.png (git commit 6cfb9e6 message mentions an MC web #3 kernel ladder on Oct 23 1:08 AM CT, earlier file not md5-matched)',17,17,[],'Red + Monte-Carlo web #3 smoothed with kernels R = 0.5..14 Gpc. Builds on #13 and #1.'),
('Harmonic decomposition 2x2',['recon_n3_2x2.png'],'Sun Oct 25 2026, 12:15 AM CT','file creation (modified again 8:23 AM CT)',25,13,[W+'cosmology/figures/fourier_resid/recon_n3_2x2.py'],
 'A and B data with model before/after removing n = 1-3 harmonics; residuals with fitted modes. Depends on #11. A chi2 25.54 to 11.46, B chi2 13.68 to 10.24.'),
]

RM='/workspace/coepdf/book/remade/';RC='/workspace/coepdf/book/recovered/'
GT='/workspace/coe-cosmology-gh/draft/current/figure_sources/src/'
SIM={'f02':(0.370,0.938),'f11':(0.779,0.944),'f12':(0.628,0.936),'f13':(0.690,0.937),'f19mc3b_stack':(0.498,0.911),'Fig06_overlay_one':(0.693,0.778),'sky_angle':(0.695,0.891),'sky_mapA_al':(0.721,0.923),'sky_mapB_al':(0.731,0.928),'sky_tab_only':(0.322,0.816)}
ST={'f06.png':('original (box)',[RC+'f06.py'],'redweb2/bodies3.py output red_web_s2000_bodies_1pct.png is md5-identical to f06.png (SSIM 1.000).'),
 'Fig01a_geometry_top_test.png':('original (git) + manual crop',[],'Per FIGURE_MANIFEST.md (coe-cosmology-gh): top/middle panels are manual crops of Fig01_geometry.png written by universe_thin.py via figs_letter_exact_python/makers/make_fig01_02_geometry_achart.py (sources in draft/current/figure_sources). Crop step not scripted.'),
 'f02.png':('remade',[RM+'f02_remade.py'],''),'f12.png':('remade',[RM+'f12_remade.py'],''),'f11.png':('remade',[RM+'f11_remade.py'],''),'f13.png':('remade',[RM+'f13_remade.py'],''),
 'sky_angle.png':('remade (spot geometry from recovered box /tmp code)',[RM+'sky_angle_remade.py',RC+'sky_mapA_al.py',RC+'spots.py'],'sky maps: recovered (box /tmp) r45c.py, axes cropped; angle plot remade.'),
 'sky_tab_only.png':('recovered (box /tmp)',[RC+'sky_tab_only.py'],'Table section of /tmp/comb45c.py; same content, 5 px wider at dpi 120.'),
 'Fig06_overlay_one.png':('remade',[RM+'Fig06_overlay_one_remade.py'],''),'f19mc3b_stack.png':('remade',[RM+'f19mc3b_stack_remade.py'],'Restack of f19mc3b.png made by redweb2/realweb/mc/mc_br_plot.py.')}
F=[(a,b,c,d,e,f,(ST[b[0]][1]+g if b[0] in ST else g),h+((' STATUS NOTE: '+ST[b[0]][2]) if b[0] in ST else '')) for (a,b,c,d,e,f,g,h) in F]
def status(files,src,scr):
    if files[0] in ST: return ST[files[0]][0]
    if not scr: return 'not found'
    return 'original (git)' if src.startswith('git commit') else 'original (box)'
def esc(s):
    for a,b in [('\\','\\textbackslash{}'),('&','\\&'),('%','\\%'),('$','\\$'),('#','\\#'),('_','\\_'),('{','\\{'),('}','\\}'),('~','\\textasciitilde{}'),('^','\\textasciicircum{}')]:
        s=s.replace(a,b) if a!='\\' else s.replace(a,'\x00')
    return s.replace('\x00','\\textbackslash{}')
TW=6.5;TH=9.0
out=[r'''\documentclass[11pt]{article}
\usepackage[letterpaper,margin=1in]{geometry}
\usepackage{fontspec,graphicx,fvextra,xcolor,amsmath}\usepackage[hidelinks]{hyperref}
\setmonofont{DejaVu Sans Mono}[Scale=0.8]
\setlength{\parindent}{0pt}\setlength{\parskip}{4pt}
\fvset{breaklines=true,breakanywhere=true,fontsize=\footnotesize,frame=single,numbers=left,numbersep=4pt}
\title{CoE paper figure book\\\large every figure, in the order it was made, with its Python source}
\author{}\date{Built Sun Oct 25 2026}
\newcommand{\flag}[1]{\textsf{[unverified: #1]}}\newcommand{\why}[1]{\textit{Validity.} #1}\usepackage{xurl}\sloppy
\begin{document}\maketitle
Order: by ``originally made'' time (CT), oldest first; ties broken by git commit, then file time to the second, then paper number. ``Copy time'' means only when the file was copied onto the box (Fri Oct 23 2026, 8:20--8:31 AM CT); the true original is likely earlier. ``This morning'' = paper as of 8:30 AM CT Oct 25; ``now'' = rewritten paper. Every image is set at \texttt{\textbackslash textwidth}; very tall images get a taller page instead of being shrunk.
\tableofcontents\clearpage''']
os.makedirs('code',exist_ok=True)
rep=[]
for i,(nm,files,dt,src,mo,now,scr,notes) in enumerate(F,1):
    out.append(r'\clearpage\section*{%d. %s}\addcontentsline{toc}{section}{%d. %s}'%(i,esc(nm),i,esc(nm)))
    out.append(r'\textbf{Originally made:} %s \quad (\textit{%s})\\ \textbf{Paper figure number:} this morning Fig.~%d; now Fig.~%d\\ \textbf{Files:} \texttt{%s}'%(esc(dt),esc(src),mo,now,esc(', '.join(files))))
    out.append(r'\\ \textbf{Code status:} \fbox{%s}'%esc(status(files,src,scr)))
    c=capmorn.get(files[0]) or capnow.get(files[0],'')
    out.append(r'\subsection*{Notes}'+esc(notes)+r'\par\textbf{Paper caption (caveats as stated):} '+re.sub(r'\\ref\{([^}]*)\}',lambda m:'{}'+esc(m.group(1)),c))
    for f in files:
        w,h=Image.open(FIG+f).size;need=TW*h/w
        g=r'\noindent\includegraphics[width=\textwidth]{%s}\par\centerline{\footnotesize\texttt{%s}}'%(FIG+f,esc(f))
        if need>TH-0.5:
            ph=need+3.0
            out.append(r'\clearpage\begingroup\pdfpageheight=%.2fin\paperheight=%.2fin\textheight=%.2fin\vsize=\textheight\makeatletter\@colht=\textheight\@colroom=\textheight\makeatother '%(ph,ph,ph-2.2)+g+r'\clearpage\endgroup ')
        else: out.append(r'\par\medskip'+g)
        b=f.rsplit('.',1)[0]
        if os.path.exists(RM+b+'_remade.png'):
            ss,mm=SIM[b]
            out.append(r'\par\medskip\noindent\textbf{REMADE (not original)} -- similarity to original: SSIM %.3f, 1-MAE %.3f\\ \includegraphics[width=0.48\textwidth,height=0.4\textheight,keepaspectratio]{%s}\hfill\includegraphics[width=0.48\textwidth,height=0.4\textheight,keepaspectratio]{%s}\par\centerline{\footnotesize left: original (unchanged) \quad right: remade}'%(ss,mm,FIG+f,RM+b+'_remade.png'))
    out.append(r'\subsection*{Generating script}')
    if not scr: out.append(r'\textbf{Generating script not found.} Searched all \texttt{.py} files under /workspace (excluding virtualenvs) for the image filename and its known aliases, and the git repos.')
    for s in scr:
        k='code/%02d_%s'%(i,os.path.basename(s));open(k,'w').write(open(s,encoding='utf-8',errors='replace').read().replace('\t','    '))
        out.append(r'\textbf{\texttt{%s}}\VerbatimInput{%s}'%(esc(s),k))
    rep.append((i,nm,status(files,src,scr),'y' if scr else ('partial' if False else 'n')))
out.append(r'\end{document}')
open('book.tex','w').write('\n'.join(out))
for r in rep:print(*r,sep=' | ')
