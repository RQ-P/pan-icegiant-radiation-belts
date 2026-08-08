import os
import sys
import json
import time
import math
import datetime
import numpy as np
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
from ice_giant_field import URANUS_G, URANUS_H, NEPTUNE_G, NEPTUNE_H
from drift_shell import mag_frame, magnetic_axis, dipole_moment, trace_line
from shell_walk import build_Lstar_table, walk_shell, from_cyl
from scipy.interpolate import PchipInterpolator
import module1c_drift_invariant as m1c
import module3_dll_star as m3
_TWO_PI = 2.0 * math.pi
CASES = [('dipole', {(1, 0): 30000.0}, {}, 1, 36, 44), ('uranus_g3', URANUS_G, URANUS_H, 3, 36, 44), ('neptune_g3', NEPTUNE_G, NEPTUNE_H, 3, 36, 44), ('neptune_o8', NEPTUNE_G, NEPTUNE_H, 8, 24, 12)]
D0, DN, DL0 = (1e-07, 3.0, 6.0)
DS = 0.02

def _np(path, **kw):
    np.savez(path, **kw)

def _meta_add(mode, entry):
    path = os.path.join(_HERE, 'meta.json')
    meta = {}
    if os.path.exists(path):
        with open(path, 'r', encoding='utf-8') as f:
            meta = json.load(f)
    meta.setdefault(mode, []).append(entry)
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(meta, f, ensure_ascii=False, indent=2)

def _table_of(name, G, H, n, n_azim, n_max_rung, ds=DS):
    frame = mag_frame(magnetic_axis(G, H))
    M1 = dipole_moment(G, H)
    tab = build_Lstar_table(G, H, n, frame, M1, ds=ds, n_azim=n_azim, n_max_rung=n_max_rung)
    return (frame, M1, tab)

def run_tables():
    for name, G, H, n, na, nr in CASES:
        t0 = time.time()
        frame, M1, tab = _table_of(name, G, H, n, na, nr)
        rungs = tab['rungs']
        Bm = np.array([r['B_m'] for r in rungs if r.get('ok')], float)
        ok = [r for r in rungs if r.get('ok')]
        Ls = np.array(tab['Ls'], float)
        L = (abs(M1) / Bm) ** (1.0 / 3.0)
        clos = np.array([r['closure'] for r in ok], float)
        asym = np.array([r['asym'] for r in ok], float)
        thN = np.array([(r['thN_deg'][0], r['thN_deg'][1]) for r in ok], float)
        _np(os.path.join(_HERE, f'lstars_{name}.npz'), B_m=Bm, Lstar=Ls, L=L, closure=clos, asym_mag=asym, thN_min_deg=thN[:, 0], thN_max_deg=thN[:, 1], M1=M1, n_max=n, n_azim=na, ds=DS, n_max_rung=nr, n_rungs=len(Bm))
        tr = m3.transform_to_Lstar(tab, M1)
        DLL = m3.D_LL(tr['L'], D0=D0, n=DN, L0=DL0)
        DLstLst = DLL * tr['J']
        dL_rel = tr['Lstar'] / tr['L'] - 1.0
        _np(os.path.join(_HERE, f'module3_{name}.npz'), B_m=tr['B_m'], L=tr['L'], Lstar=tr['Lstar'], dLstar_dL=tr['dLstar_dL'], J=tr['J'], dL_rel=dL_rel, DLL=DLL, DLstLst=DLstLst, fold_k=tr['fold_k'], main_k=tr['main_k'], M1=M1, D0=D0, n=DN, L0=DL0)
        _meta_add('tables', dict(case=name, n_max=n, n_azim=na, n_max_rung=nr, M1=float(M1), n_rungs=int(len(Bm)), wall_s=round(time.time() - t0, 1)))

def run_crosscheck():
    for name, G, H, n, B_ms in [('uranus_g3', URANUS_G, URANUS_H, 3, [120.0, 400.0, 1200.0]), ('neptune_g3', NEPTUNE_G, NEPTUNE_H, 3, [120.0, 400.0, 800.0])]:
        for na, ds in [(36, 0.02), (72, 0.01)]:
            t0 = time.time()
            rc = m1c.crosscheck_case(name, G, H, n, B_ms, n_azim=na, ds=ds)
            rows = np.array(rc['rows'], float)
            ds_tag = f'{int(ds * 1000):03d}'
            _np(os.path.join(_HERE, f'crosscheck_{name}_a{na}_d{ds_tag}.npz'), B_m=rows[:, 0], L_orbit=rows[:, 1], L_sw=rows[:, 2], L_naive=rows[:, 3], rel=rows[:, 4], n_good=rows[:, 5].astype(int), M1=rc['M1'], n_azim=na, ds=ds)
            _meta_add('crosscheck', dict(case=name, n_azim=na, ds=ds, wall_s=round(time.time() - t0, 1)))

def run_conv():
    plans = {'uranus_g3': [(36, 0.02), (72, 0.02), (72, 0.01)], 'neptune_g3': [(36, 0.02), (72, 0.02)]}
    bm_grid = [120.0, 400.0, 800.0, 1200.0]
    for name, G, H, n in [('uranus_g3', URANUS_G, URANUS_H, 3), ('neptune_g3', NEPTUNE_G, NEPTUNE_H, 3)]:
        if name not in _planet_filter():
            continue
        res_rows = []
        n_azims, ds_list = ([], [])
        for na, ds in plans[name]:
            t0 = time.time()
            _, M1, tab = _table_of(name, G, H, n, na, 44, ds=ds)
            tr = m3.transform_to_Lstar(tab, M1)
            clos = [r['closure'] for r in tab['rungs'] if r.get('ok')]
            L_grid = np.array([float(tab['f_L'](math.log(b))) for b in bm_grid])
            n_azims.append(na)
            ds_list.append(ds)
            res_rows.append(dict(n_rungs=len(tab['Bs']), Bs=np.asarray(tab['Bs'], float), Ls=np.asarray(tab['Ls'], float), L_grid=L_grid, fold_Bm=np.asarray(tr['B_m'][tr['fold_k']], float), fold_Lstar=np.asarray(tr['Lstar'][tr['fold_k']], float), closure_max=float(max(clos)) if clos else float('nan'), wall_s=round(time.time() - t0, 1)))
        _np(os.path.join(_HERE, f'conv_{name}.npz'), n_azims=np.array(n_azims), ds=np.array(ds_list), n_rungs=np.array([r['n_rungs'] for r in res_rows]), Bs=np.array([r['Bs'] for r in res_rows], dtype=object), Ls=np.array([r['Ls'] for r in res_rows], dtype=object), L_grid=np.array([r['L_grid'] for r in res_rows]), fold_Bm=np.array([r['fold_Bm'] for r in res_rows], dtype=object), fold_Lstar=np.array([r['fold_Lstar'] for r in res_rows], dtype=object), closure_max=np.array([r['closure_max'] for r in res_rows]))
        _meta_add('conv', dict(case=name, bm_grid=bm_grid))

def run_fascan():
    planets = _planet_filter()
    alphas = _alpha_filter()
    limits = {('neptune_g3', 45.0): dict(max_steps=250000, sample_every=3), ('neptune_g3', 30.0): dict(max_steps=150000, sample_every=3)}
    for name, G, H, n in [('uranus_g3', URANUS_G, URANUS_H, 3), ('neptune_g3', NEPTUNE_G, NEPTUNE_H, 3)]:
        if name not in planets:
            continue
        t0 = time.time()
        frame, M1, tab = _table_of(name, G, H, n, 36, 44)
        al, L90, Lfin, split, relA, covs = ([], [], [], [], [], [])
        dth, ovr, spr, lim = ([], [], [], [])
        for a in alphas:
            lim_conf = limits.get((name, a), {})
            r = m1c.finite_alpha_phi_conservation(name, G, H, n, 400.0, alpha_eq_deg=a, q_eff=4.0, revs=3, max_steps=lim_conf.get('max_steps', 600000), sample_every=lim_conf.get('sample_every', 1), tab=tab, verbose=False)
            al.append(a)
            L90.append(r['L90'])
            Lfin.append(r['L_finite'])
            split.append(r['split'])
            relA.append(r['relA'])
            covs.append(r['cov'])
            dth.append(r['dth_max'])
            ovr.append(r['overshoot_rel'])
            spr.append(r['b_eq_spread'])
            lim.append(r['limited'])
        _np(os.path.join(_HERE, f'fascan_{name}.npz'), alpha=np.array(al), L90=np.array(L90), L_finite=np.array(Lfin), split=np.array(split), relA=np.array(relA), cov=np.array(covs), dth_max_deg=np.array(dth), overshoot_rel=np.array(ovr), b_eq_spread=np.array(spr), limited=np.array(lim, bool), sample_every=np.array([limits.get((name, a), {}).get('sample_every', 1) for a in al]), max_steps=np.array([limits.get((name, a), {}).get('max_steps', 600000) for a in al]), B_m_seed=400.0, M1=M1, q_eff=4.0, revs=3)
        _meta_add('fascan', dict(case=name, B_m_seed=400.0, wall_s=round(time.time() - t0, 1)))

def run_topology():
    targets = {'uranus_g3': [4.0, 5.5, 7.0], 'neptune_g3': [3.0, 4.5, 6.0]}
    n_line = 12
    for name, G, H, n in [('uranus_g3', URANUS_G, URANUS_H, 3), ('neptune_g3', NEPTUNE_G, NEPTUNE_H, 3)]:
        t0 = time.time()
        frame, M1, tab = _table_of(name, G, H, n, 36, 44)
        fL = tab['f_L']
        Bg = np.logspace(np.log10(tab['B_lo']), np.log10(tab['B_hi']), 400)
        shells = []
        for Ltar in targets[name]:
            i = int(np.argmin(np.abs(fL(np.log(Bg)) - Ltar)))
            B_sol = float(Bg[i])
            L_sol = float(fL(math.log(B_sol)))
            r = m1c._find_rho_eq(B_sol, 0.0, G, H, n, frame)
            if r is None:
                r = m1c._find_rho_eq(B_sol, 0.3, G, H, n, frame)
            if r is None:
                continue
            seed = (r, 0.0, 0.0)
            sh = walk_shell(B_sol, seed, G, H, n, frame, ds=DS, n_azim=36)
            if not sh['ok']:
                continue
            lines = []
            idx = np.linspace(0, len(sh['pts']) - 1, n_line).astype(int)
            for j in idx:
                rho_c, zeta_c, phi_c = sh['pts'][j]
                pos = from_cyl(rho_c, zeta_c, phi_c, frame)
                t = trace_line(pos, G, H, n, ds=DS)
                if t['ok'] and t['path'] is not None:
                    lines.append(np.asarray(t['path'], float))
            shells.append(dict(L_target=Ltar, B_m=B_sol, L_sol=L_sol, thN=np.asarray(sh['thN']), phN=np.asarray(sh['phN']), thS=np.asarray(sh['thS']), phS=np.asarray(sh['phS']), pts=np.asarray([from_cyl(*p, frame) for p in sh['pts']]), closure=sh['closure'], lines=lines))
        _np(os.path.join(_HERE, f'topology_{name}.npz'), axis=np.asarray(frame, float), shells=np.array(shells, dtype=object), M1=M1)
        _meta_add('topology', dict(case=name, n_shells=len(shells), wall_s=round(time.time() - t0, 1)))

def _branch_coeffs(case):
    z = np.load(os.path.join(_HERE, f'module3_{case}.npz'))
    mk = z['main_k']
    L = z['L'][mk]
    Ls = z['Lstar'][mk]
    J = z['J'][mk]
    o = np.argsort(Ls)
    L, Ls, J = (L[o], Ls[o], J[o])
    D0, n0, L0 = (float(z['D0']), float(z['n']), float(z['L0']))
    Linv = PchipInterpolator(Ls, L)
    DLL = D0 * (Linv(Ls) / L0) ** n0
    Dcorr = DLL * J
    Dnaive = D0 * (Ls / L0) ** n0
    return dict(Ls=Ls, L=L, J=J, Linv=Linv, Dcorr=Dcorr, Dnaive=Dnaive, M1=float(z['M1']), D0=D0, n=n0, L0=L0)

def _fp_mean_time(Ls, D, x_tar, x_a):
    from scipy.integrate import solve_bvp
    sel = Ls <= x_tar
    xs = Ls[sel]
    Ds = np.asarray(D)[sel]
    lnD = PchipInterpolator(xs, np.log(Ds))
    Dp = PchipInterpolator(xs, Ds)

    def ode(x, y):
        return np.vstack([y[1], -lnD.derivative()(x) * y[1] - 1.0 / Dp(x)])

    def bc(ya, yb):
        return np.array([ya[1], yb[0]])
    y0 = np.zeros((2, xs.size))
    y0[0] = (xs - x_tar) ** 2
    y0[1] = 2.0 * (xs - x_tar)
    sol = solve_bvp(ode, bc, xs, y0, max_nodes=6000)
    if not sol.success:
        raise RuntimeError(f'FP mean-time solve_bvp failed: {sol.message}')
    tau = np.full(Ls.size, np.nan)
    tau[sel] = sol.sol(xs)[0]
    return tau

def _mc_first_pass(Ls, D, x_inj, x_tar, n_p=10000, n_step=400000, seed=7):
    rng = np.random.default_rng(seed)
    Dp = PchipInterpolator(Ls, D)
    dD = Dp.derivative()
    x_a, x_b = (float(Ls[0]), float(Ls[-1]))
    Dmax = float(D.max())
    dt = 20.0 * (x_b - x_a) ** 2 / (n_step * Dmax)
    x = np.full(n_p, x_inj)
    tau = np.full(n_p, np.nan)
    done = np.zeros(n_p, bool)
    for k in range(n_step):
        x = x + dD(x) * dt + np.sqrt(2.0 * Dp(x) * dt) * rng.standard_normal(n_p)
        x = np.minimum(x, x_b)
        x = np.maximum(x, x_a)
        arr = ~done & (x >= x_tar)
        if arr.any():
            tau[arr] = k * dt
            done[arr] = True
        if done.all():
            break
    return tau

def run_mcsim():
    for case, Linj, Ltar in [('uranus_g3', 3.0, 6.9), ('neptune_g3', 2.9, 7.0)]:
        t0 = time.time()
        S = _branch_coeffs(case)
        Ls = S['Ls']
        tau_c = _mc_first_pass(Ls, S['Dcorr'], Linj, Ltar)
        tau_n = _mc_first_pass(Ls, S['Dnaive'], Linj, Ltar)
        fp_c = _fp_mean_time(Ls, S['Dcorr'], Ltar, Ls[0])
        fp_n = _fp_mean_time(Ls, S['Dnaive'], Ltar, Ls[0])
        i_inj = int(np.argmin(np.abs(Ls - Linj)))
        fold_vals = np.load(os.path.join(_HERE, f'module3_{case}.npz'))['Lstar']
        fold_inside = np.asarray([float(f) for f in np.load(os.path.join(_HERE, f'module3_{case}.npz'))['Lstar'][np.load(os.path.join(_HERE, f'module3_{case}.npz'))['fold_k']]])
        fold_inside = fold_inside[(fold_inside > Linj) & (fold_inside < Ltar)]
        mc_c = float(np.nanmean(tau_c))
        mc_n = float(np.nanmean(tau_n))
        r = mc_c / max(fp_c[i_inj], 1e-300)
        rn = mc_n / max(fp_n[i_inj], 1e-300)
        _np(os.path.join(_HERE, f'mcsim_{case}.npz'), Ls=Ls, L=S['L'], J=S['J'], Dcorr=S['Dcorr'], Dnaive=S['Dnaive'], tau_corr=tau_c, tau_naive=tau_n, fp_corr=fp_c, fp_naive=fp_n, Linj=Linj, Ltar=Ltar, fold_inside=fold_inside, mc_mean_corr=mc_c, mc_mean_naive=mc_n, fp_mean_corr=float(fp_c[i_inj]), fp_mean_naive=float(fp_n[i_inj]))
        _meta_add('mcsim', dict(case=case, Linj=Linj, Ltar=Ltar, wall_s=round(time.time() - t0, 1)))

def _planet_filter():
    try:
        i = sys.argv.index('--planet')
        return set(sys.argv[i + 1].split(','))
    except (ValueError, IndexError):
        return {'uranus_g3', 'neptune_g3'}

def _alpha_filter():
    try:
        i = sys.argv.index('--alpha')
        return [float(a) for a in sys.argv[i + 1].split(',')]
    except (ValueError, IndexError):
        return [75.0, 60.0, 45.0, 30.0]

def main():
    modes = [a for a in sys.argv[1:] if a.startswith('--') and a not in ('--planet', '--alpha')]
    if not modes or '--all' in modes:
        modes = ['--tables', '--crosscheck', '--conv', '--fascan', '--topology', '--mcsim']
    for m in modes:
        t0 = time.time()
        {'--tables': run_tables, '--crosscheck': run_crosscheck, '--conv': run_conv, '--fascan': run_fascan, '--topology': run_topology, '--mcsim': run_mcsim}[m]()
if __name__ == '__main__':
    main()