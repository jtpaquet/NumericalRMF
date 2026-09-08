"""Shared matplotlib style, imported as `from utils import plot_options`."""
from cycler import cycler

plot_options = {
    "axes.grid" : False,
    'axes.prop_cycle' : cycler('color', 'krbgcmy'),
    'axes.unicode_minus' : False, # Otherwise, bug with negative number in axes
    'errorbar.capsize' : 3,
    'figure.dpi':100,
    'figure.figsize':(12,9),
    'font.family' : 'serif',
    'font.sans-serif': 'Computer Modern Sans Serif',
    'font.serif': 'cmr10',
    'font.size': 20,
    'image.cmap' : 'jet',
    'lines.color' : 'k',
    'lines.linewidth' : 0.7,
    'lines.markerfacecolor': 'w',
    'lines.markersize': 10,
    'markers.fillstyle': 'none',
    'mathtext.fontset':'cm',
    'scatter.edgecolors' : 'face',
    'scatter.marker' : 'o',
    'text.usetex' : False,
    'xtick.direction' : 'in',
    'xtick.major.size' : 4.0,
    'xtick.minor.visible' : True,
    'xtick.top': True,
    'ytick.direction' : 'in',
    'ytick.major.size' : 4.0,
    'ytick.minor.visible' : True,
    'ytick.right': True,
    }
