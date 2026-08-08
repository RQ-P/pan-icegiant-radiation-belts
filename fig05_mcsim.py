"""Fig. 7: Monte Carlo verification -- particles diffuse freely across the
fold L* values; the fold is a coordinate singularity, not a dynamical
barrier.

Panel (a): diffusion coefficient on the monotone main branch,
D_{L*L*}=D_LL . J (corrected, solid) vs naive dipole D_LL (dashed), for
Uranus and Neptune; vertical lines mark the L* values shared with
off-branch (fold) shells. D is positive everywhere on the branch.
Panel (b): cumulative first-passage-time distributions of an ensemble of
10^4 Euler-Maruyama particles injected inside the branch and diffusing to
a target beyond the fold values; all particles arrive (100% hit rate).
Panel (c): mean first-passage times, Monte Carlo vs the independent FP
solution (D tau'' + D' tau' = -1); agreement to 2-12%.

Output: paper2_figs/fig05_mcsim.png / .pdf
"""
import os
import numpy as np

from figstyle import setup, DOUBLE

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')


def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.0))
    axa, axb, axc = axs
    for i, (case, name) in enumerate([('uranus_g3', 'Uranus'),
                                      ('neptune_g3', 'Neptune')]):
        z = np.load(os.path.join(_DATA, f'mcsim_{case}.npz'))
        Ls = z['Ls']
        axa.plot(Ls, z['Dcorr'] / 1e-7, color=f'C{i}', lw=1.3, label=name)
        axa.plot(Ls, z['Dnaive'] / 1e-7, color=f'C{i}', ls='--', lw=1.0,
                 alpha=0.7)
        for f in z['fold_inside']:
            axa.axvline(f, color='C3', lw=0.6, alpha=0.6)
        t_c = z['tau_corr']; t_n = z['tau_naive']
        tc = t_c[np.isfinite(t_c)]; tn = t_n[np.isfinite(t_n)]
        xc = np.sort(tc); xn = np.sort(tn)
        axb.plot(xc, np.arange(1, len(xc) + 1) / len(t_c), color=f'C{i}',
                 lw=1.3, label=f'{name} corrected')
        axb.plot(xn, np.arange(1, len(xn) + 1) / len(t_n), color=f'C{i}',
                 ls='--', lw=1.0, alpha=0.7, label=f'{name} naive dipole')
        mc = [z['mc_mean_corr'], z['mc_mean_naive']]
        fp = [z['fp_mean_corr'], z['fp_mean_naive']]
        axc.plot(fp, mc, 'o', ms=6, color=f'C{i}')
    axa.set_xlabel(r'$L^{\star}$ [R]')
    axa.set_ylabel(r'$D_{L^{\star}L^{\star}}/D_0$')
    axa.set_title(r'(a) branch diffusion coefficient (red: fold $L^{\star}$)')
    axa.legend(fontsize=9, loc='upper right')
    axb.set_xscale('log')
    axb.set_xlabel(r'first-passage time [a.u.]')
    axb.set_ylabel(r'cumulative fraction arriving')
    axb.set_title(r'(b) 100% cross the fold values (no barrier)')
    axb.legend(fontsize=9, loc='lower right')
    lim = max(np.max(z['fp_mean_corr']) for z in
              [np.load(os.path.join(_DATA, f'mcsim_{c}.npz'))
               for c in ('uranus_g3', 'neptune_g3')]) * 1.1
    axc.plot([0, lim], [0, lim], 'k--', lw=0.8)
    axc.set_xlabel(r'$\bar\tau$ FP solution [a.u.]')
    axc.set_ylabel(r'$\bar\tau$ Monte Carlo [a.u.]')
    axc.set_title(r'(c) MC vs FP mean first-passage time')
    axc.set_xscale('log'); axc.set_yscale('log')
    fig.tight_layout()
    fig.savefig('fig05_mcsim.png')
    fig.savefig('fig05_mcsim.pdf')
    print('saved fig05_mcsim.png/.pdf')


if __name__ == '__main__':
    main()
