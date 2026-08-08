import os
import sys
import numpy as np
from scipy.interpolate import PchipInterpolator
_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
if _ROOT not in sys.path:
    sys.path.insert(0, _ROOT)
import module3_dll_star as m3
D0, DN, DL0 = (1e-07, 3.0, 6.0)

def main():
    for case in ['uranus_g3', 'neptune_g3']:
        conv = np.load(os.path.join(_HERE, f'conv_{case}.npz'), allow_pickle=True)
        old = np.load(os.path.join(_HERE, f'module3_{case}.npz'))
        i = next((k for k in range(len(conv['n_azims'])) if conv['n_azims'][k] == 72 and abs(conv['ds'][k] - 0.02) < 1e-09))
        Bs = np.asarray(conv['Bs'][i], float)
        Ls = np.asarray(conv['Ls'][i], float)
        tab = dict(f_L=PchipInterpolator(np.log(Bs), Ls), Bs=Bs, Ls=Ls)
        M1 = float(old['M1'])
        tr = m3.transform_to_Lstar(tab, M1)
        DLL = m3.D_LL(tr['L'], D0=D0, n=DN, L0=DL0)
        DLstLst = DLL * tr['J']
        dL_rel = tr['Lstar'] / tr['L'] - 1.0
        np.savez(os.path.join(_HERE, f'module3_{case}_a36.npz'), **{k: old[k] for k in old.files})
        np.savez(os.path.join(_HERE, f'module3_{case}.npz'), B_m=tr['B_m'], L=tr['L'], Lstar=tr['Lstar'], dLstar_dL=tr['dLstar_dL'], J=tr['J'], dL_rel=dL_rel, DLL=DLL, DLstLst=DLstLst, fold_k=tr['fold_k'], main_k=tr['main_k'], M1=M1, D0=D0, n=DN, L0=DL0, n_azim=72, ds=0.02)
        Jm = tr['J'][tr['main_k']]
if __name__ == '__main__':
    main()