import os
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.integrate import solve_bvp
from figstyle import setup, DOUBLE
_HERE = os.path.dirname(os.path.abspath(__file__))
_DATA = os.path.join(_HERE, '..', 'paper2_data')
D0, N, L0 = (1e-07, 3.0, 6.0)
KAPPA2 = 0.5

def d_coeffs(case):
    z = np.load(os.path.join(_DATA, f'module3_{case}.npz'))
    mk = z['main_k']
    L, Ls, J = (z['L'][mk], z['Lstar'][mk], z['J'][mk])
    o = np.argsort(Ls)
    L, Ls, J = (L[o], Ls[o], J[o])
    Linv = PchipInterpolator(Ls, L)
    Jp = PchipInterpolator(Ls, J)
    Lstar = np.linspace(Ls.min(), Ls.max(), 200)
    D_naive = D0 * (Lstar / L0) ** N
    D_corr = D0 * (Linv(Lstar) / L0) ** N * Jp(Lstar)
    return (Lstar, D_naive, D_corr, Jp, Ls)

def solve_f(Lstar, Dfun):
    Dp = PchipInterpolator(Lstar, Dfun)
    dlnD = PchipInterpolator(Lstar, np.log(Dfun))
    fa, fb = (1.0, 0.1)

    def ode(x, y):
        return np.vstack([y[1], -dlnD.derivative()(x) * y[1] + KAPPA2 * D0 * y[0] / Dp(x)])

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
    axa.plot(Lstar, D_corr / D0, 'C0', lw=1.3, label='corrected: $D_{LL}\\cdot J$')
    axa.set_xlabel('$L^{\\star}$ [R]')
    axa.set_ylabel('$D(L^{\\star})/D_0$')
    axa.set_yscale('log')
    axa.set_title('(a) Uranus $D_{L^{\\star}L^{\\star}}$ remapping')
    axa.legend(fontsize=7)
    axb.plot(Lstar, f_naive, 'k--', lw=1.0, label='naive dipole')
    axb.plot(Lstar, f_corr, 'C0', lw=1.3, label='corrected ($J$ remap)')
    axb.set_xlabel('$L^{\\star}$ [R]')
    axb.set_ylabel('$f(L^{\\star})$ [normalized]')
    axb.set_title('(b) steady-state $f$ (Paper I formalism)')
    axb.legend(fontsize=7)
    fig.tight_layout()
    fig.savefig('fig07_bridge_transport.png')
    fig.savefig('fig07_bridge_transport.pdf')
if __name__ == '__main__':
    main()