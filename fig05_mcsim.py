import os
import numpy as np
from figstyle import setup, DOUBLE
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')

def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.0))
    axa, axb, axc = axs
    for i, (case, name) in enumerate([('uranus_g3', 'Uranus'), ('neptune_g3', 'Neptune')]):
        z = np.load(os.path.join(_DATA, f'mcsim_{case}.npz'))
        Ls = z['Ls']
        axa.plot(Ls, z['Dcorr'] / 1e-07, color=f'C{i}', lw=1.3, label=name)
        axa.plot(Ls, z['Dnaive'] / 1e-07, color=f'C{i}', ls='--', lw=1.0, alpha=0.7)
        for f in z['fold_inside']:
            axa.axvline(f, color='C3', lw=0.6, alpha=0.6)
        t_c = z['tau_corr']
        t_n = z['tau_naive']
        tc = t_c[np.isfinite(t_c)]
        tn = t_n[np.isfinite(t_n)]
        xc = np.sort(tc)
        xn = np.sort(tn)
        axb.plot(xc, np.arange(1, len(xc) + 1) / len(t_c), color=f'C{i}', lw=1.3, label=f'{name} corrected')
        axb.plot(xn, np.arange(1, len(xn) + 1) / len(t_n), color=f'C{i}', ls='--', lw=1.0, alpha=0.7, label=f'{name} naive dipole')
        mc = [z['mc_mean_corr'], z['mc_mean_naive']]
        fp = [z['fp_mean_corr'], z['fp_mean_naive']]
        axc.plot(fp, mc, 'o', ms=6, color=f'C{i}')
    axa.set_xlabel('$L^{\\star}$ [R]')
    axa.set_ylabel('$D_{L^{\\star}L^{\\star}}/D_0$')
    axa.set_title('(a) branch diffusion coefficient (red: fold $L^{\\star}$)')
    axa.legend(fontsize=6, loc='upper right')
    axb.set_xscale('log')
    axb.set_xlabel('first-passage time [a.u.]')
    axb.set_ylabel('cumulative fraction arriving')
    axb.set_title('(b) 100% cross the fold values (no barrier)')
    axb.legend(fontsize=6, loc='lower right')
    lim = max((np.max(z['fp_mean_corr']) for z in [np.load(os.path.join(_DATA, f'mcsim_{c}.npz')) for c in ('uranus_g3', 'neptune_g3')])) * 1.1
    axc.plot([0, lim], [0, lim], 'k--', lw=0.8)
    axc.set_xlabel('$\\bar\\tau$ FP solution [a.u.]')
    axc.set_ylabel('$\\bar\\tau$ Monte Carlo [a.u.]')
    axc.set_title('(c) MC vs FP mean first-passage time')
    axc.set_xscale('log')
    axc.set_yscale('log')
    fig.tight_layout()
    fig.savefig('fig05_mcsim.png')
    fig.savefig('fig05_mcsim.pdf')
if __name__ == '__main__':
    main()