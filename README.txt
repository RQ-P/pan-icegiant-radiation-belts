================================================================
Drift-shell folding and radial transport modification in Uranus
and Neptune magnetospheres caused by realistic non-dipolar
magnetic fields (Icarus)
Reproducibility package - README
================================================================

This package reproduces all data tables and figures of the paper.

1. ENVIRONMENT
----------------------------------------------------------------
- Python 3.13 with numpy, scipy, matplotlib.
- Figures use the local SciencePlots style
  (SciencePlots/SciencePlots-master); paper2_figs/figstyle.py loads the
  "science" style with the "no-latex" override (no TeX required).

2. LAYERS (computation and plotting are decoupled)
----------------------------------------------------------------
- paper2_data/build_data.py : THE compute layer. Runs the frozen core
  (ice_giant_field.py / drift_shell.py / shell_walk.py, read-only
  imports) plus module1c_drift_invariant.py and module3_dll_star.py,
  and stores every result as .npz under paper2_data/.
      python paper2_data/build_data.py --tables      (L* tables + J)
      python paper2_data/build_data.py --crosscheck  (Path A, two res.)
      python paper2_data/build_data.py --conv        (resolution scan)
      python paper2_data/build_data.py --fascan      (splitting vs alpha)
      python paper2_data/build_data.py --topology    (field-line point
                                                      clouds for fig 1)
      python paper2_data/build_data.py --mcsim       (Monte Carlo
                                                      crossing test, fig 7)
  meta.json records case parameters and wall times.

- paper2_figs/fig0X_*.py     : THE plotting layer. Read only the .npz
  files; no physics is recomputed. Each writes PNG and PDF.

3. DATA FILES
----------------------------------------------------------------
  lstars_<case>.npz          B_m, Lstar, L, closure, asym_mag,
                             north-footpoint colatitude range
  module3_<case>.npz         L, Lstar, dLstar_dL, J, dL_rel, DLL,
                             DLstLst, fold_k, main_k
  crosscheck_<case>_a<na>_d<ds>.npz
                             orbit vs shell-walk vs naive rows
  conv_<case>.npz            full tables + L*(B_m grid) at 2-3
                             resolutions (fold convergence)
  fascan_<case>.npz            split vs alpha + reliability diagnostics
  mcsim_<case>.npz             Monte Carlo first-passage across fold
                               values (Euler-Maruyama on the main branch,
                               FP cross-check): tau_corr/tau_naive,
                               fp_corr/fp_naive, hit rates, fold_inside
  topology_<case>.npz          traced field lines + drift-shell point
                               clouds for Fig. 1(a,b)
  meta.json                  run parameters and wall times

4. SCRIPT-TO-FIGURE CORRESPONDENCE
----------------------------------------------------------------
  Figure | Script                       | Verified claim
  -------+-----------------------------+--------------------------------
  Fig. 1 | fig01_crosscheck.py          | orbit vs shell-walk agreement
         |                              | (0.5-5.8%); fold-affected point
  Fig. 2 | fig02_lstar_tables.py        | dipole L-shell failure
         |                              | (L*/L-1 up to -45%)
  Fig. 3 | fig03_split_vs_alpha.py      | drift-shell splitting vs alpha;
         |                              | loss-cone-limited points hollow
  Fig. 4 | fig04_topology.py            | computed field topology (a,b);
         |                              | L*(L) folding, J=0 (c)
  Fig. 5 | fig05_mcsim.py               | Monte Carlo: 100% of particles
         |                              | cross the fold L* values; mean
         |                              | first-passage time matches the
         |                              | FP solution to 2-12% (folds are
         |                              | coordinate singularities, not
         |                              | dynamical barriers)
  Fig. 6 | fig06_jacobian_folds.py      | L*(L) folding (a), dL*/dL
         |                              | crossing zero (b), remapped D with
         |                              | dipole reference (c), resolution
         |                              | convergence + fold counts (d)
  Fig. 7 | fig07_bridge_transport.py    | steady-state f(L*) naive vs
         |                              | J-corrected on the main branch

5. KEY NUMBERS (as stated in the paper; converged n_azim=72, ds=0.02
   tables unless noted)
----------------------------------------------------------------
- Uranus deg-3: L* in [1.86, 7.73] R, L*/L-1 in [-45.3%, -8.5%],
  main-branch J in [0.130, 2.354], 5 fold points
  (B_m 14.7/17.9/39.7/48.5/1424 nT).
- Neptune deg-3: L* in [2.28, 10.18] R, L*/L-1 in [-26.3%, -1.3%],
  main-branch J in [0.041, 1.818], 4 fold points (incl. loss-cone
  boundary at B_m ~ 780 nT, n_azim=72).
- module3_<case>.npz hold the converged tables; module3_<case>_a36.npz
  are the coarser n_azim=36 backups.
- Splitting at B_m=400 nT: Uranus alpha=45 deg +12.2% (reliable);
  Neptune alpha=75 deg +0.8% (reliable), lower alphas loss-cone-limited.
- Monte Carlo (fig 7): 100% of particles cross the fold L* values
  (Uranus fold values 5.71/6.73/6.87 inside the branch; Neptune 3.13);
  MC mean first-passage time matches the FP solution to 2-12%; the
  coordinate factor dilates the time scale by ~1.25 (Uranus) / ~1.0
  (Neptune) -- folds are coordinate singularities, not dynamical
  barriers.
- Dipole self-checks: V1/V2 ~1e-8, V4 ~1e-15, finite-alpha Phi
  conservation ~2e-6.



