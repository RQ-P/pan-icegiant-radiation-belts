"""Fig. 3: Cross-validation of the drift-shell reconstruction (Path A).

Panels (a,b): L*_orbit (independent |B|=B_m equatorial orbit), L*_shellwalk
(constant-B_min shell walking) and naive dipole L for the Uranus (a) and
Neptune (b) B_m grid, at n_azim=36/ds=0.02 and n_azim=72/ds=0.01.
Panels (c,d): relative orbit-vs-shellwalk deviation at both resolutions;
the B_m=800 nT point for Neptune sits inside the fold region and is marked.

Output: paper2_figs/fig01_crosscheck.png / .pdf
"""
import os
import numpy as np

from figstyle import setup, DOUBLE

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')


def load(case, na, ds):
    tag = f'a{na}_d{int(ds * 1000):03d}'
    return np.load(os.path.join(_DATA, f'crosscheck_{case}_{tag}.npz'))


def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(2, 2, figsize=(7.5, 6.2))
    for j, (case, name) in enumerate([('uranus_g3', 'Uranus (deg-3)'),
                                      ('neptune_g3', 'Neptune (deg-3)')]):
        axv = axs[0, j]
        axr = axs[1, j]
        for (na, ds), ms, lbl in [((36, 0.02), 6, r'$n_{\rm az}$=36'),
                                  ((72, 0.01), 8, r'$n_{\rm az}$=72')]:
            z = load(case, na, ds)
            B = z['B_m']
            axv.plot(B, z['L_orbit'], 'o', ms=ms, label=lbl + ' orbit')
            axv.plot(B, z['L_sw'], 's', ms=ms, mfc='none', label=lbl + ' shell-walk')
            axv.plot(B, z['L_naive'], '^', ms=ms, mfc='none', alpha=0.7,
                     label=lbl + ' naive dipole')
            axr.plot(B, z['rel'] * 100.0, 'o-', ms=ms, label=lbl)
        axv.set_xscale('log')
        axv.set_yscale('log')
        axv.set_xlabel(r'$B_m$ [nT]')
        axv.set_ylabel(r'$L^{\star}$ [R]')
        axv.set_title(f'({chr(97 + j)}) {name}: L* methods')
        axv.legend(fontsize=9, )
        axr.axhline(0, color='k', lw=0.6)
        axr.set_xscale('log')
        axr.set_xlabel(r'$B_m$ [nT]')
        axr.set_ylabel(r'$(L^{\star}_{\rm orbit}/L^{\star}_{\rm sw}-1)$ [%]')
        if case == 'neptune_g3':
            axr.annotate('fold region\n(B_m=800)', xy=(800, -19),
                         xytext=(420, -40), fontsize=9, 
                         arrowprops=dict(arrowstyle='->', lw=0.7))
        axr.set_title(f'({chr(99 + j)}) {name}: orbit vs shell-walk')
        axr.legend(fontsize=9, )
    fig.tight_layout()
    fig.savefig('fig01_crosscheck.png')
    fig.savefig('fig01_crosscheck.pdf')
    print('saved fig01_crosscheck.png/.pdf')


if __name__ == '__main__':
    main()
