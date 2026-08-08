import os
import numpy as np
from figstyle import setup, DOUBLE
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')

def main():
    plt = setup(DOUBLE)
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for i, (case, name) in enumerate([('uranus_g3', 'Uranus (deg-3)'), ('neptune_g3', 'Neptune (deg-3)')]):
        z = np.load(os.path.join(_DATA, f'fascan_{case}.npz'))
        a, s = (z['alpha'], z['split'] * 100.0)
        lim = z['limited']
        ax.plot(a[~lim], s[~lim], 'o-', ms=5, color=f'C{i}', label=name + ' (reliable)')
        ax.plot(a[lim], s[lim], 'o--', ms=5, mfc='none', color=f'C{i}', label=name + ' (loss-cone limited)')
    ax.axhline(0, color='k', lw=0.6)
    ax.set_xlabel('$\\alpha_{\\rm eq,seed}$ [deg]')
    ax.set_ylabel('drift-shell splitting $\\Delta L^{\\star}/L^{\\star}$ [%]')
    ax.set_title('Finite pitch angle: $L^{\\star}(\\alpha)/L^{\\star}_{90^\\circ}-1$')
    ax.legend(fontsize=7, loc='best')
    ax.set_xticks([30, 45, 60, 75])
    fig.tight_layout()
    fig.savefig('fig03_split_vs_alpha.png')
    fig.savefig('fig03_split_vs_alpha.pdf')
if __name__ == '__main__':
    main()