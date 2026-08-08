"""Fig. 6: Bridge to Paper I — steady-state radial diffusion in the L*
coordinate (Module 4).

The Paper I formalism (steady-state 1D radial diffusion with loss,
f'' + (D'/D) f' - f/(D tau) = 0) is solved as a two-point BVP on the
monotone main branch of Uranus. The naive dipole treatment uses
D_naive(L*) = D0 (L*/L0)^n; the corrected treatment remaps
D(L*) = D_LL(L(L*)) . (dL*/dL)^2 with L(L*) the inverse of the main-branch
mapping, which locally suppresses transport where the Jacobian dips
(loss-cone/fold-adjacent region).

Panel (a): D(L*) naive vs corrected (log).
Panel (b): steady-state f(L*) (normalized), naive vs corrected.

Output: paper2_figs/fig07_bridge_transport.png / .pdf
"""
import os
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_bvp

from figstyle import setup, DOUBLE

_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')

D0, N, L0 = 1.0e-7, 3.0, 6.0
KAPPA2 = 0.5          # 1/(D0 tau), Paper I convention; illustrative


def d_coeffs(case):
    z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
    mk = z['main_k']
    L, Ls, J = z['L'][mk], z['Lstar'][mk], z['J'][mk]
    o = np.argsort(Ls)
    L, Ls, J = L[o], Ls[o], J[o]
    Linv = PchipInterpolator(Ls, L)                    # L(L*), main branch
    Jp = PchipInterpolator(Ls, J)                      # J(L*)
    Lstar = np.linspace(Ls.min(), Ls.max(), 200)
    D_naive = D0 * (Lstar / L0) ** N
    D_corr = D0 * (Linv(Lstar) / L0) ** N * Jp(Lstar)
    return Lstar, D_naive, D_corr, Jp, Ls


def solve_f(Lstar, Dfun):
    Dp = PchipInterpolator(Lstar, Dfun)
    dlnD = PchipInterpolator(Lstar, np.log(Dfun))
    fa, fb = 1.0, 0.1

    def ode(x, y):
        # f'' + (D'/D) f' - kappa^2 D0 f / D = 0,  kappa^2 = 1/(D0 tau)
        return np.vstack([y[1], -dlnD.derivative()(x) * y[1]
                          + KAPPA2 * D0 * y[0] / Dp(x)])

    def bc(ya, yb):
        return np.array([ya[0] - fa, yb[0] - fb])

    y0 = np.zeros((2, Lstar.size))
    y0[0] = fa + (fb - fa) * (Lstar - Lstar[0]) / (Lstar[-1] - Lstar[0])
    y0[1] = (fb - fa) / (Lstar[-1] - Lstar[0])
    sol = solve_bvp(ode, bc, Lstar, y0, max_nodes=5000)
    if not sol.success:
        raise RuntimeError(f'solve_bvp failed: {sol.message}')
    return sol.sol(Lstar)[0]


def main():
    plt = setup(DOUBLE)
    fig, axs = plt.subplots(1, 2, figsize=(9.5, 4.0))
    axa, axb = axs
    Lstar, D_naive, D_corr, Jp, Ls_branch = d_coeffs('uranus_g3')
    f_naive = solve_f(Lstar, D_naive)
    f_corr = solve_f(Lstar, D_corr)
    axa.plot(Lstar, D_naive / D0, 'k--', lw=1.0, label='naive dipole')
    axa.plot(Lstar, D_corr / D0, 'C0', lw=1.3,
             label=r'corrected: $D_{LL}\cdot J$')
    axa.set_xlabel(r'$L^{\star}$ [R]')
    axa.set_ylabel(r'$D(L^{\star})/D_0$')
    axa.set_yscale('log')
    axa.set_title(r'(a) Uranus $D_{L^{\star}L^{\star}}$ remapping')
    axa.legend(fontsize=7)
    axb.plot(Lstar, f_naive, 'k--', lw=1.0, label='naive dipole')
    axb.plot(Lstar, f_corr, 'C0', lw=1.3,
             label=r'corrected ($J$ remap)')
    axb.set_xlabel(r'$L^{\star}$ [R]')
    axb.set_ylabel(r'$f(L^{\star})$ [normalized]')
    axb.set_title(r'(b) steady-state $f$ (Paper I formalism)')
    axb.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig('fig07_bridge_transport.png')
    fig.savefig('fig07_bridge_transport.pdf')
    print('saved fig07_bridge_transport.png/.pdf')


if __name__ == '__main__':
    main()
