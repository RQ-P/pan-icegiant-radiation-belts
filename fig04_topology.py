"""Fig. 1: Drift-shell topology and L*(L) mapping in the ice giants.

Panel (a): realistic computed field geometry (Connerney harmonic model,
deg-3) in the magnetic-axis meridian plane (rho, zeta): traced field lines
(light), drift-shell min-B loci (colored closed curves), magnetic axis
(dashed). Data: paper2_data/topology_*.npz.
Panel (b): L*(L) mapping from the Module 3 Jacobian data; fold points
(J=0, L* coordinate non-invertible) marked with x. Data:
paper2_data/module3_*.npz.
Panel (c): concept flow diagram.

Output: paper2_figs/fig04_topology.png / .pdf
"""
import os
import numpy as np

from figstyle import setup, DOUBLE

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')


def _cyl(points):
    p = np.asarray(points, float)
    rho = np.hypot(p[..., 0], p[..., 1])
    zeta = p[..., 2]
    return rho, zeta


def panel_topology(ax, case, title):
    z = np.load(os.path.join(_DATA, f'topology_{case}.npz'),
                allow_pickle=True)
    ax.add_patch(__import__('matplotlib.patches', fromlist=['Circle']).Circle(
        (0, 0), 1.0, facecolor='#d8d8d8', edgecolor='k', lw=0.8, zorder=2))
    ax.text(0, 0, 'planet', ha='center', va='center', fontsize=9, zorder=3)
    for shell in z['shells']:
        Ls = shell['L_sol']
        for line in shell['lines']:
            r, s = _cyl(line)
            ax.plot(r, s, lw=0.35, color='#9ec1e3', alpha=0.7, zorder=1)
        pts = shell['pts']
        r, s = _cyl(pts)
        ax.plot(r, s, lw=1.4, zorder=4,
                label=fr'$L^{{\star}}\approx{Ls:.1f}$')
        ax.plot(r[0], s[0], 'o', ms=3, zorder=5)
    ax.axvline(0, color='k', ls='--', lw=0.7, zorder=1)
    ax.text(0.02, 4.6, 'magnetic axis', rotation=90, fontsize=9, va='top')
    ax.set_xlabel(r'$\rho$ [R]')
    ax.set_ylabel(r'$\zeta$ [R]')
    ax.set_xlim(-6, 6)
    ax.set_ylim(-6, 6)
    ax.set_aspect('equal')
    ax.set_title(title)
    ax.legend(fontsize=9, loc='lower left')


def panel_mapping(ax):
    zr = np.load(os.path.join(_DATA, 'module3_uranus_g3.npz'))
    zn = np.load(os.path.join(_DATA, 'module3_neptune_g3.npz'))
    for z, name, c in [(zr, 'Uranus', 0), (zn, 'Neptune', 1)]:
        L, Ls = z['L'], z['Lstar']
        mk = z['main_k']
        ax.plot(L[mk], Ls[mk], color=f'C{c}', lw=1.3,
                label=f'{name} (monotone branch)')
        nk = np.setdiff1d(np.arange(len(L)), mk)
        if nk.size:
            ax.plot(L[nk], Ls[nk], 'o', ms=2.5, mfc='none',
                    color=f'C{c}', alpha=0.8)
        if z['fold_k'].size:
            fk = z['fold_k']
            ax.plot(L[fk], Ls[fk], 'x', color='C3', ms=5, zorder=5)
    Lm = np.linspace(1, 14, 2)
    ax.plot(Lm, Lm, 'k--', lw=0.9, label='dipole reference')
    ax.set_xlabel(r'$L=(M_1/B_{\rm eq})^{1/3}$ [R]')
    ax.set_ylabel(r'$L^{\star}$ [R]')
    ax.set_title('L*(L) mapping: folds (x) mark $J=0$')
    ax.legend(fontsize=9, loc='upper left')


def panel_concept(ax):
    import matplotlib.patches as mpatches
    boxes = ['realistic non-dipolar\nharmonics (Voyager)',
             r'$L^{\star}$ folding\n(non-monotone shells)',
             r'modified diffusion\n$D_{L^{\star}L^{\star}}$',
             'localized\ntransport change']
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 3)
    ax.axis('off')
    w, h = 2.35, 1.5
    x0, y0 = 0.1, 0.75
    for i, b in enumerate(boxes):
        x = x0 + i * 2.45
        ax.add_patch(mpatches.FancyBboxPatch(
            (x, y0), w, h, boxstyle='round,pad=0.08',
            fc='#eef3f8', ec='C7', lw=1.0))
        ax.text(x + w / 2, y0 + h / 2, b, ha='center', va='center',
                fontsize=7)
        if i < len(boxes) - 1:
            ax.annotate('', xy=(x + w + 0.18, y0 + h / 2),
                        xytext=(x + w + 0.02, y0 + h / 2),
                        arrowprops=dict(arrowstyle='->', lw=1.2))
    ax.set_title('interpretation')


def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(1, 3, figsize=(11.5, 4.2))
    panel_topology(axs[0], 'uranus_g3', '(a) Uranus: computed topology')
    panel_topology(axs[1], 'neptune_g3', '(b) Neptune: computed topology')
    panel_mapping(axs[2])
    fig.tight_layout()
    fig.savefig('fig04_topology.png')
    fig.savefig('fig04_topology.pdf')
    print('saved fig04_topology.png/.pdf')


if __name__ == '__main__':
    main()
