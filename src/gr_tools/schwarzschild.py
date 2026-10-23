"""Schwarzschild solution: symbolic derivation, metric constructor and orbit helpers.

Conventions: signature (-, +, +, +), geometric units G = c = 1, coordinates (t, r, theta, phi).
"""
import itertools
import math

import sympy as sp

from .engine import GREngine

t, r, th, ph = sp.symbols("t r theta phi", real=True)
COORDS = [t, r, th, ph]


def schwarzschild_metric(M=1):
    """Return the Schwarzschild metric as a diagonal ``sp.Matrix``.

    Args:
        M: Mass parameter (a number or a SymPy symbol), in geometric units.
    """
    f = 1 - 2 * M / r
    return sp.diag(-f, 1 / f, r**2, r**2 * sp.sin(th) ** 2)


def derive_schwarzschild_metric(verbose=True):
    """Derive the Schwarzschild metric from the vacuum Einstein equations.

    Starts from the static, spherically symmetric ansatz
    ``ds^2 = -e^{2 alpha(r)} dt^2 + e^{2 beta(r)} dr^2 + r^2 dOmega^2``,
    solves ``R_mu_nu = 0`` and validates the result.

    Returns:
        dict with keys ``metric`` (sp.Matrix), ``is_vacuum`` (bool) and
        ``kretschmann`` (sp.Expr).
    """

    def say(*args):
        if verbose:
            print(*args)

    alpha, beta = sp.Function("alpha")(r), sp.Function("beta")(r)

    g_ansatz = sp.diag(-sp.exp(2 * alpha), sp.exp(2 * beta), r**2, r**2 * sp.sin(th) ** 2)
    engine = GREngine(COORDS, g_ansatz)
    Ric = engine.ricci()

    say("--- Step 1: Solving R^t_t - R^r_r = 0 ---")
    eq_alpha_beta = sp.simplify(engine.g_inv[0, 0] * Ric[0, 0] - engine.g_inv[1, 1] * Ric[1, 1])
    say(f"Relationship eq: {eq_alpha_beta} = 0")
    say("=> alpha'(r) = -beta'(r) => beta(r) = -alpha(r) (setting integration constant to 0 for asymptotic flatness)\n")

    say("--- Step 2: Solving R_thth = 0 to find alpha(r) ---")
    y = sp.Function("y")(r)

    R_thth_sub = Ric[2, 2].subs(beta, -alpha).subs(alpha, sp.log(y) / 2).doit()
    expr_y = sp.simplify(R_thth_sub)
    say(f"ODE for y(r) = exp(2*alpha): {expr_y} = 0")

    sol = sp.dsolve(sp.Eq(expr_y, 0), y)
    say(f"General solution for y(r): {sol.rhs}")

    constants = list(sol.rhs.free_symbols - {r})
    if not constants:
        raise ValueError("SymPy did not return a constant of integration.")
    rs = sp.symbols("r_s", real=True)
    y_final = sol.rhs.subs(constants[0], -rs)
    say(f"Applying physical boundary conditions (constant = -r_s): y(r) = {y_final}\n")

    say("--- Step 3: Constructing and Validating the Derived Metric ---")
    alpha_final = sp.log(y_final) / 2
    beta_final = -alpha_final

    g_derived = sp.simplify(
        sp.diag(-sp.exp(2 * alpha_final), sp.exp(2 * beta_final), r**2, r**2 * sp.sin(th) ** 2)
    )

    val_engine = GREngine(COORDS, g_derived)
    ricci_final = val_engine.ricci()
    is_vacuum = all(
        sp.simplify(ricci_final[i, j]) == 0 for i, j in itertools.product(range(4), repeat=2)
    )
    say(f"Einstein Vacuum Equations (R_mu_nu = 0) satisfied dynamically: {is_vacuum}")

    k_scalar = val_engine.kretschmann()
    if verbose:
        dt, dr, dth, dph = sp.symbols("dt dr dtheta dphi")
        ds2 = (
            g_derived[0, 0] * dt**2
            + g_derived[1, 1] * dr**2
            + g_derived[2, 2] * dth**2
            + g_derived[3, 3] * dph**2
        )
        print("\nLine Element (ds^2):")
        sp.pprint(sp.Eq(sp.Symbol("ds^2"), ds2))
        print("\nKretschmann Scalar (K = R_abcd R^abcd):")
        sp.pprint(k_scalar)

    return {"metric": g_derived, "is_vacuum": is_vacuum, "kretschmann": k_scalar}


def effective_potential(radius, L, M=1.0):
    """Effective potential per unit mass for timelike geodesics.

    ``(dr/dtau)^2 = E^2 - V_eff(r)`` with ``V_eff = (1 - 2M/r)(1 + L^2/r^2)``.
    """
    return (1 - 2 * M / radius) * (1 + L**2 / radius**2)


def circular_orbit_angular_momentum(radius, M=1.0):
    """Specific angular momentum of a timelike circular orbit (needs ``radius > 3M``)."""
    return (M * radius**2 / (radius - 3 * M)) ** 0.5


def circular_orbit_initial_conditions(radius, M=1.0):
    """Equatorial circular-orbit initial state ``(x, u)`` with ``u^mu u_mu = -1``."""
    ut = (1 - 3 * M / radius) ** -0.5
    uphi = (M / radius**3) ** 0.5 * ut
    return [0.0, radius, math.pi / 2, 0.0], [ut, 0.0, 0.0, uphi]


def conserved_quantities(x, u, M=1.0):
    """Energy ``E = -u_t`` and angular momentum ``L = u_phi`` per unit mass."""
    _, radius, theta, _ = x
    E = (1 - 2 * M / radius) * u[0]
    L = radius**2 * math.sin(theta) ** 2 * u[3]
    return float(E), float(L)


if __name__ == "__main__":
    derive_schwarzschild_metric()
