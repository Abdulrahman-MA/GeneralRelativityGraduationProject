"""Killing vector test: ``L_xi g = 0``."""
import sympy as sp

from .engine import GREngine


def lie_derivative_of_metric(engine: GREngine, xi):
    """Return the matrix ``(L_xi g)_{mu nu}`` for the vector field ``xi^mu``.

    ``(L_xi g)_{mu nu} = xi^l d_l g_{mu nu} + g_{l nu} d_mu xi^l + g_{mu l} d_nu xi^l``.
    """
    n, x, g = engine.n, engine.coords, engine.g
    out = sp.zeros(n, n)
    for m in range(n):
        for v in range(n):
            out[m, v] = sp.simplify(
                sum(xi[l] * sp.diff(g[m, v], x[l]) for l in range(n))
                + sum(g[l, v] * sp.diff(xi[l], x[m]) for l in range(n))
                + sum(g[m, l] * sp.diff(xi[l], x[v]) for l in range(n))
            )
    return out


def is_killing(engine: GREngine, xi):
    """True if ``xi`` (list of components in the engine's coordinates) is a Killing vector."""
    return lie_derivative_of_metric(engine, xi) == sp.zeros(engine.n, engine.n)
