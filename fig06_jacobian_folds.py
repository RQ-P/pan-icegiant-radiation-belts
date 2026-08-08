import os
import numpy as np
from figstyle import setup, DOUBLE
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')
IG = [('uranus_g3', 'Uranus (deg-3)'), ('neptune_g3', 'Neptune (deg-3)')]

def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(2, 2, figsize=(9.5, 7.6))
    (axa, axb), (axc, axd) = axs
    zd = np.load(os.path.join(_DATA, 'module3_dipole.npz'))
    Lm = np.linspace(1, 14, 2)
    axa.plot(Lm, Lm, 'k--', lw=0.9, label='dipole reference')
    for i, (case, name) in enumerate(IG):
        z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
        mk = z['main_k']
        axa.plot(z['L'][mk], z['Lstar'][mk], lw=1.3, label=name)
        nk = np.setdiff1d(np.arange(len(z['L'])), mk)
        if nk.size:
            axa.plot(z['L'][nk], z['Lstar'][nk], 'o', ms=2.5, mfc='none', color=f'C{i}', alpha=0.8)
        if z['fold_k'].size:
            fk = z['fold_k']
            axa.plot(z['L'][fk], z['Lstar'][fk], 'x', color='C3', ms=6, zorder=5)
    axa.set_xlabel('$L=(M_1/B_{\\rm eq})^{1/3}$ [R]')
    axa.set_ylabel('$L^{\\star}$ [R]')
    axa.set_title('(a) $L^{\\star}(L)$ folds; x: $\\mathcal{J}=0$')
    axa.legend(fontsize=6, loc='upper left')
    axb.axhline(0, color='C3', lw=1.6)
    axb.text(1.05, 0.02, 'coordinate fold  ($\\mathcal{J}=0$)', fontsize=6, color='C3')
    for i, (case, name) in enumerate(IG):
        z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
        mk = z['main_k']
        d = z['dLstar_dL']
        axb.plot(z['Lstar'][mk], d[mk], lw=1.3, label=name)
        nk = np.setdiff1d(np.arange(len(z['Lstar'])), mk)
        if nk.size:
            axb.plot(z['Lstar'][nk], d[nk], 'o', ms=2.5, mfc='none', color=f'C{i}', alpha=0.8)
        if z['fold_k'].size:
            for k in z['fold_k']:
                axb.axvline(z['Lstar'][k], color='C7', lw=0.6, alpha=0.5)
    axb.set_xlabel('$L^{\\star}$ [R]')
    axb.set_ylabel('$\\mathrm{d}L^{\\star}/\\mathrm{d}L$')
    axb.set_title('(b) derivative crosses zero at the folds')
    axb.legend(fontsize=6, loc='upper right')
    bm_grid = [120.0, 400.0, 800.0, 1200.0]
    for i, (case, name) in enumerate(IG):
        z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
        mk = z['main_k']
        Ls, D = (z['Lstar'], z['DLstLst'])
        D0, n0, L0 = (float(z['D0']), float(z['n']), float(z['L0']))
        Ddip = D0 * (Ls / L0) ** n0
        axc.plot(Ls, Ddip, ls='--', lw=0.9, color='C6', alpha=0.8)
        ln, = axc.plot(Ls[mk], D[mk], lw=1.3, label=name)
        nk = np.setdiff1d(np.arange(len(Ls)), mk)
        if nk.size:
            axc.plot(Ls[nk], D[nk], 'o', ms=2.5, ls='none', color=ln.get_color(), alpha=0.8)
        if z['fold_k'].size:
            for k in z['fold_k']:
                axc.axvline(z['Lstar'][k], color='C7', lw=0.6, alpha=0.5)
    axc.plot([], [], ls='--', color='C6', lw=0.9, label='dipole $D_{LL}$')
    ymax = max((np.load(os.path.join(_DATA, f'module3_{c}.npz'))['DLstLst'].max() for c, _ in IG))
    axc.set_ylim(0, ymax * 1.05)
    axc.set_xlabel('$L^{\\star}$ [R]')
    axc.set_ylabel('$D_{L^{\\star}L^{\\star}}$ [a.u.]')
    axc.set_title('(c) remapped $D$: suppression at fold-adjacent $L^{\\star}$')
    axc.legend(fontsize=6, loc='lower right')
    for i, (case, name) in enumerate(IG):
        z = np.load(os.path.join(_DATA, f'conv_{case}.npz'), allow_pickle=True)
        labels = [f'{int(a)}/{ds:.2f}' for a, ds in zip(z['n_azims'], z['ds'])]
        x = np.arange(len(labels))
        lg = z['L_grid']
        for j, b in enumerate(bm_grid):
            axd.plot(x, lg[:, j], 'o-', ms=4, lw=1.0, color=f'C{i}', alpha=0.45 + 0.55 * (1 - j / len(bm_grid)))
        folds = [len(np.asarray(f)) for f in z['fold_Bm']]
        axd.annotate(f'{name}\nfolds: {folds}', xy=(x[0], lg[0, 1]), xytext=(x[0] + 0.1, lg[0, 1] + 1.5), fontsize=6, va='center')
    axd.set_xticks(range(len(labels)))
    axd.set_xticklabels(labels, rotation=15, fontsize=6)
    axd.set_xlabel('resolution (n_azim / ds)')
    axd.set_ylabel('$L^{\\star}(B_m)$ [R]')
    axd.set_title('(d) table convergence; fold counts per resolution')
    axd.set_ylim(0, 7)
    fig.tight_layout()
    fig.savefig('fig06_jacobian_folds.png')
    fig.savefig('fig06_jacobian_folds.pdf')
if __name__ == '__main__':
    main()