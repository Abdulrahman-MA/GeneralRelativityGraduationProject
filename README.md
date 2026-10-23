# General Relativity: From Flat Spacetime to the Schwarzschild Solution

Bachelor of Science graduation project (supervised by Prof. Dr. Elsayed Lashin). It summarises and expands the proofs in the first five chapters of Sean Carroll's *Spacetime and Geometry*, then adds symbolic and numerical Python tools that check the main results.

> The thesis follows Carroll's book. It is an independent study aid, not a replacement for the book.

## What's in this repository

| Folder | Contents |
|--------|----------|
| [`thesis/`](thesis) | The final thesis (49 pages): `thesis.tex` and the compiled [`thesis.pdf`](thesis/thesis.pdf) |
| [`chapter_notes/`](chapter_notes) | The full-length notes the thesis was condensed from, one folder per chapter, each with `.tex` and `.pdf` |
| [`src/gr_tools/`](src/gr_tools) | A small Python package for symbolic curvature, Killing vectors and geodesics |
| [`tests/`](tests) | Tests that check the code against known results (see below) |
| [`notebooks/`](notebooks) | Derivations and plots used in the chapters |
| [`scripts/`](scripts) | 3D visualisations of Killing vector fields (PyVista) |
| [`presentation/`](presentation) | Slides and Birkhoff-theorem animations |
| [`assets/figures/`](assets/figures) | Figures used in the notes and notebooks |

## Thesis outline

1. Special Relativity and Flat Spacetime
2. Manifolds
3. Curvature
4. Gravitation (Einstein's equation, the Palatini formulation)
5. The Schwarzschild Solution (Birkhoff's theorem, coordinate vs. curvature singularities)

## Chapter notes

| Chapter | Folder | PDF |
|---------|--------|-----|
| 1. Special Relativity | [`ch1_special_relativity`](chapter_notes/ch1_special_relativity) | [`ch1.pdf`](chapter_notes/ch1_special_relativity/ch1.pdf) |
| 2. Manifolds | [`ch2_manifolds`](chapter_notes/ch2_manifolds) | [`ch2.pdf`](chapter_notes/ch2_manifolds/ch2.pdf) |
| 3. Curvature | [`ch3_curvature`](chapter_notes/ch3_curvature) | [`ch3.pdf`](chapter_notes/ch3_curvature/ch3.pdf) |
| 4. Gravitation | [`ch4_gravitation`](chapter_notes/ch4_gravitation) | [`ch4.pdf`](chapter_notes/ch4_gravitation/ch4.pdf) |
| 5. The Schwarzschild Solution | [`ch5_schwarzschild`](chapter_notes/ch5_schwarzschild) | [`ch5.pdf`](chapter_notes/ch5_schwarzschild/ch5.pdf) |
| 6. Causality (partial: orientability) | [`ch6_causality`](chapter_notes/ch6_causality) | [`ch6.pdf`](chapter_notes/ch6_causality/ch6.pdf) |

## Python tools (`gr_tools`)

```bash
python -m venv .venv
.venv\Scripts\activate            # Linux/macOS: source .venv/bin/activate
pip install -e ".[dev]"           # add ",viz" for the PyVista scripts, ",notebooks" for JupyterLab
pytest
```

```python
import sympy as sp
from gr_tools import GREngine, schwarzschild_metric, derive_schwarzschild_metric
from gr_tools.schwarzschild import COORDS

M = sp.symbols("M", positive=True)
engine = GREngine(COORDS, schwarzschild_metric(M))
engine.ricci()          # zero matrix: vacuum solution
engine.kretschmann()    # 48*M**2/r**6

derive_schwarzschild_metric()   # derives the metric from the vacuum Einstein equations, step by step
```

| Module | What it does |
|--------|--------------|
| `engine.py` | `GREngine`: Christoffel symbols, Riemann, Ricci and Kretschmann scalar for any metric |
| `schwarzschild.py` | Symbolic derivation of the Schwarzschild metric, effective potential, circular-orbit helpers |
| `killing.py` | Tests whether a vector field is a Killing vector via the Lie derivative of the metric |
| `geodesics.py` | Numerical geodesic integration (DOP853) using the engine's Christoffel symbols |

### What the tests check

- A flat plane in polar coordinates has zero curvature, and the unit 2-sphere has `R_mn = g_mn`.
- The Schwarzschild metric is Ricci-flat, with Kretschmann scalar `48 M^2 / r^6`.
- The derivation from the ansatz recovers `-(1 - r_s/r)` and `K = 12 r_s^2 / r^6`.
- The three rotation fields of the 2-sphere and `d/dt`, `d/dphi` of Schwarzschild are Killing vectors. `d/dtheta` and `d/dr` are not.
- Along numerically integrated geodesics, `u.u = -1`, the energy `E` and the angular momentum `L` stay constant, and a circular orbit at `r = 10M` stays circular.
- The effective potential is extremal at the circular-orbit angular momentum, and `L^2 = 12 M^2` at the ISCO.

## Rebuilding the PDFs

The sources use TikZ, and the chapter notes use `\includesvg`, which needs `-shell-escape` and Inkscape.

```bash
cd thesis && pdflatex thesis.tex && pdflatex thesis.tex
cd chapter_notes/ch1_special_relativity && pdflatex -shell-escape ch1.tex
```

## Author

Abdelrahman Mohamed Anwar
