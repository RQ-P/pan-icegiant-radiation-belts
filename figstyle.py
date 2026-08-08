import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
_HERE = os.path.dirname(os.path.abspath(__file__))
_STYLES = os.path.normpath(os.path.join(_HERE, '..', 'SciencePlots', 'SciencePlots-master', 'src', 'scienceplots', 'styles'))
_BASE = [os.path.join(_STYLES, 'science.mplstyle'), os.path.join(_STYLES, 'misc', 'no-latex.mplstyle')]
SINGLE = (3.37, 3.37)
DOUBLE = (6.69, 4.6)
CYCLE = ['#0C5DA5', '#00B945', '#FF9500', '#FF2C00', '#845B97', '#474747', '#9e9e9e']

def setup(figsize=DOUBLE, dpi=600, rows=1):
    if figsize is SINGLE:
        figsize = (3.37, 3.0)
    elif figsize is DOUBLE:
        figsize = (6.69, 4.6)
    plt.style.use(_BASE)
    plt.rcParams.update({'figure.figsize': figsize, 'savefig.dpi': dpi, 'figure.dpi': 100, 'font.size': 9, 'axes.labelsize': 9, 'axes.titlesize': 9, 'xtick.labelsize': 8, 'ytick.labelsize': 8, 'legend.fontsize': 8, 'lines.linewidth': 1.2, 'axes.grid': False, 'savefig.bbox': 'tight', 'savefig.pad_inches': 0.03, 'axes.prop_cycle': matplotlib.cycler(color=CYCLE)})
    return plt