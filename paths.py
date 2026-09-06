"""Output locations, anchored to the repo root so scripts run from anywhere.

    results/    pickles and tabulated data
    figures/    pdf and png
    figures/animations/   mp4
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
ANIMATIONS = FIGURES / 'animations'


def results(name):
    """Absolute path for a file under results/, creating the directory."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    return str(RESULTS / name)


def figure(name):
    """Absolute path for a file under figures/, creating the directory."""
    p = FIGURES / name
    p.parent.mkdir(parents=True, exist_ok=True)
    return str(p)


def animation(name):
    """Absolute path for a file under figures/animations/."""
    ANIMATIONS.mkdir(parents=True, exist_ok=True)
    return str(ANIMATIONS / name)
