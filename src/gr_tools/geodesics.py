"""Numerical geodesic integration for any metric handled by :class:`GREngine`."""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

from .engine import GREngine


def four_velocity_norm(engine: GREngine, x, u):
    """Return ``g_{mu nu} u^mu u^nu`` at the point ``x`` (numeric metric needed)."""
    g = np.array(sp.lambdify(engine.coords, engine.g, "numpy")(*x), dtype=float)
    u = np.asarray(u, dtype=float)
    return float(u @ g @ u)


def integrate_geodesic(engine: GREngine, x0, u0, lam_span, num=None, rtol=1e-10, atol=1e-12):
    """Integrate ``d^2x^mu/dlam^2 = -Gamma^mu_{nu rho} u^nu u^rho``.

    Args:
        engine: engine built from a metric without free symbols (e.g. ``M`` already a number).
        x0: initial coordinates.
        u0: initial four-velocity components ``dx^mu/dlam``.
        lam_span: ``(lam_start, lam_end)`` affine parameter interval.
        num: if given, number of equally spaced output samples.
        rtol, atol: integrator tolerances (DOP853).

    Returns:
        The ``scipy`` result object; ``sol.y[:n]`` are the coordinates and
        ``sol.y[n:]`` the four-velocity along the curve.
    """
    n = engine.n
    gamma = sp.lambdify(engine.coords, engine.christoffel().tolist(), "numpy")

    def rhs(_, y):
        x, u = y[:n], y[n:]
        G = np.array(gamma(*x), dtype=float)
        return np.concatenate([u, -np.einsum("mnr,n,r->m", G, u, u)])

    t_eval = None if num is None else np.linspace(*lam_span, num)
    y0 = np.concatenate([np.asarray(x0, dtype=float), np.asarray(u0, dtype=float)])
    return solve_ivp(rhs, lam_span, y0, method="DOP853", t_eval=t_eval, rtol=rtol, atol=atol)
