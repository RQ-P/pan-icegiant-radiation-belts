import os
import numpy as np
from figstyle import setup, DOUBLE
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')
CASES = [('dipole', 'Dipole (ref)'), ('uranus_g3', 'Uranus (deg-3)'), ('neptune_g3', 'Neptune (deg-3)'), ('neptune_o8', 'Neptune (O8)')]

def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.0))
    axa, axb = axs
    for i, (case, name) in enumerate(CASES):
        if case == 'dipole':
            z = np.load(os.path.join(_DATA, f'lstars_{case}.npz'))
            axa.loglog(z['B_m'], z['Lstar'], 'o-', ms=3, lw=1.0, label=name)
            continue
        z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
        axa.loglog(z['B_m'], z['Lstar'], 'o-', ms=3, lw=1.0, label=name)
        axa.loglog(z['B_m'], z['L'], ls='--', lw=0.9, alpha=0.7, color=f'C{i}')
    axa.set_xlabel('$B_m$ [nT]')
    axa.set_ylabel('$L^{\\star}$ [R]')
    axa.set_title('(a) $L^{\\star}(B_m)$ tables; dashed: dipole reference $L$')
    axa.legend(fontsize=7, loc='lower left')
    for i, (case, name) in enumerate(CASES):
        if case == 'dipole':
            continue
        z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
        Ls, mk = (z['Lstar'], z['main_k'])
        d = z['dL_rel']
        axb.plot(Ls[mk], d[mk] * 100.0, lw=1.3, label=name)
        nk = np.setdiff1d(np.arange(len(Ls)), mk)
        if nk.size:
            axb.plot(Ls[nk], d[nk] * 100.0, 'o', ms=2.5, mfc='none', color=f'C{i}', alpha=0.8)
        if z['fold_k'].size:
            fk = z['fold_k']
            axb.plot(Ls[fk], d[fk] * 100.0, 'x', color='C3', ms=5, zorder=5)
    axb.axhline(0, color='k', lw=0.6)
    axb.set_xlabel('$L^{\\star}$ [R]')
    axb.set_ylabel('$L^{\\star}/L-1$ [%]')
    axb.set_title('(b) dipole-L failure; x: fold points ($J=0$)')
    axb.legend(fontsize=7, loc='lower right')
    fig.tight_layout()
    fig.savefig('fig02_lstar_tables.png')
    fig.savefig('fig02_lstar_tables.pdf')
if __name__ == '__main__':
    main()