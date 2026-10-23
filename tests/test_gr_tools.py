import math

import numpy as np
import pytest
import sympy as sp

from gr_tools import (
    GREngine,
    circular_orbit_angular_momentum,
    circular_orbit_initial_conditions,
    conserved_quantities,
    derive_schwarzschild_metric,
    effective_potential,
    four_velocity_norm,
    integrate_geodesic,
    is_killing,
    schwarzschild_metric,
)
from gr_tools.schwarzschild import COORDS, r, ph, t, th

M = sp.symbols("M", positive=True)


@pytest.fixture(scope="module")
def schwarzschild_engine():
    return GREngine(COORDS, schwarzschild_metric(M))


@pytest.fixture(scope="module")
def sphere_engine():
    return GREngine([th, ph], sp.diag(1, sp.sin(th) ** 2))


@pytest.fixture(scope="module")
def numeric_schwarzschild():
    return GREngine(COORDS, schwarzschild_metric(1))


def test_polar_plane_is_flat():
    engine = GREngine([r, th], sp.diag(1, r**2))
    gamma = engine.christoffel()
    assert sp.simplify(gamma[0, 1, 1] + r) == 0  # Gamma^r_{theta theta} = -r
    assert sp.simplify(gamma[1, 0, 1] - 1 / r) == 0  # Gamma^theta_{r theta} = 1/r
    assert all(sp.simplify(x) == 0 for x in sp.flatten(engine.riemann().tolist()))


def test_unit_sphere_ricci_equals_metric(sphere_engine):
    assert sp.simplify(sphere_engine.ricci() - sphere_engine.g) == sp.zeros(2, 2)


def test_schwarzschild_is_ricci_flat(schwarzschild_engine):
    assert schwarzschild_engine.ricci() == sp.zeros(4, 4)


def test_schwarzschild_kretschmann(schwarzschild_engine):
    assert sp.simplify(schwarzschild_engine.kretschmann() - 48 * M**2 / r**6) == 0


def test_derivation_recovers_schwarzschild():
    result = derive_schwarzschild_metric(verbose=False)
    rs = sp.symbols("r_s", real=True)
    g = result["metric"]
    assert result["is_vacuum"]
    assert sp.simplify(g[0, 0] + (1 - rs / r)) == 0
    assert sp.simplify(g[1, 1] - 1 / (1 - rs / r)) == 0
    assert sp.simplify(result["kretschmann"] - 12 * rs**2 / r**6) == 0


def test_killing_vectors_of_unit_sphere(sphere_engine):
    cot = sp.cos(th) / sp.sin(th)
    fields = [[0, 1], [-sp.sin(ph), -cot * sp.cos(ph)], [sp.cos(ph), -cot * sp.sin(ph)]]
    assert all(is_killing(sphere_engine, xi) for xi in fields)
    assert not is_killing(sphere_engine, [1, 0])  # d/dtheta is not Killing


def test_killing_vectors_of_schwarzschild(schwarzschild_engine):
    assert is_killing(schwarzschild_engine, [1, 0, 0, 0])  # d/dt
    assert is_killing(schwarzschild_engine, [0, 0, 0, 1])  # d/dphi
    assert not is_killing(schwarzschild_engine, [0, 1, 0, 0])  # d/dr


def test_effective_potential_extremum_at_circular_orbit():
    radius, h = 10.0, 1e-5
    L = circular_orbit_angular_momentum(radius)
    slope = (effective_potential(radius + h, L) - effective_potential(radius - h, L)) / (2 * h)
    assert abs(slope) < 1e-8


def test_isco_angular_momentum():
    assert circular_orbit_angular_momentum(6.0) ** 2 == pytest.approx(12.0)


def test_circular_orbit_stays_circular(numeric_schwarzschild):
    x0, u0 = circular_orbit_initial_conditions(10.0)
    assert four_velocity_norm(numeric_schwarzschild, x0, u0) == pytest.approx(-1.0)
    sol = integrate_geodesic(numeric_schwarzschild, x0, u0, (0, 500), num=200)
    assert sol.success
    assert np.allclose(sol.y[1], 10.0, atol=1e-6)
    assert sol.y[3, -1] > 10  # phi keeps advancing


def test_eccentric_orbit_conserves_invariants(numeric_schwarzschild):
    radius, L = 20.0, 4.2
    f = 1 - 2 / radius
    uphi = L / radius**2
    ut = math.sqrt((1 + radius**2 * uphi**2) / f)
    x0, u0 = [0.0, radius, math.pi / 2, 0.0], [ut, 0.0, 0.0, uphi]

    sol = integrate_geodesic(numeric_schwarzschild, x0, u0, (0, 3000), num=300)
    assert sol.success
    E0, L0 = conserved_quantities(x0, u0)
    for k in range(sol.y.shape[1]):
        x, u = sol.y[:4, k], sol.y[4:, k]
        assert four_velocity_norm(numeric_schwarzschild, x, u) == pytest.approx(-1.0, abs=1e-7)
        E, Lk = conserved_quantities(x, u)
        assert E == pytest.approx(E0, abs=1e-7)
        assert Lk == pytest.approx(L0, abs=1e-7)
    assert sol.y[1].min() < radius  # the orbit actually moves radially
