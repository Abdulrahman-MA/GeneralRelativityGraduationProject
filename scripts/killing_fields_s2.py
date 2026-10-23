"""Plot the three Killing vector fields of the unit 2-sphere (the angular part of Schwarzschild).

Usage:  python scripts/killing_fields_s2.py {z,x,y} [--out figure.png]

In coordinates (theta, phi) the rotation generators are
    z-rotation:  K = d/dphi
    x-rotation:  K = -sin(phi) d/dtheta - cot(theta) cos(phi) d/dphi
    y-rotation:  K = +cos(phi) d/dtheta - cot(theta) sin(phi) d/dphi
and ``gr_tools.is_killing`` verifies them (see tests/test_gr_tools.py).
Requires the ``viz`` extra:  pip install -e ".[viz]"
"""
import argparse

import numpy as np
import pyvista as pv

# (K^theta, K^phi) as functions of (theta, phi), in the coordinate basis.
FIELDS = {
    "z": lambda th, ph: (np.zeros_like(th), np.ones_like(th)),
    "x": lambda th, ph: (-np.sin(ph), -np.cos(ph) / np.tan(th)),
    "y": lambda th, ph: (np.cos(ph), -np.sin(ph) / np.tan(th)),
}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("axis", choices=FIELDS)
    parser.add_argument("--out", default=None, help="screenshot path (default: show a window)")
    args = parser.parse_args()

    theta, phi = np.meshgrid(np.linspace(0.01, np.pi - 0.01, 20), np.linspace(0, 2 * np.pi, 60))
    x, y, z = np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)

    e_theta = np.stack([np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)], axis=-1)
    e_phi = np.stack([-np.sin(phi), np.cos(phi), np.zeros_like(phi)], axis=-1)

    k_theta, k_phi = FIELDS[args.axis](theta, phi)
    # Coordinate basis -> orthonormal basis: |d/dphi| = sin(theta).
    V = k_theta[..., None] * e_theta + (np.sin(theta) * k_phi)[..., None] * e_phi
    V /= np.linalg.norm(V, axis=-1, keepdims=True)

    points = np.c_[x.ravel(), y.ravel(), z.ravel()]
    cloud = pv.PolyData(points)
    cloud["vectors"] = V.reshape(-1, 3)
    arrows = cloud.glyph(orient="vectors", scale=False, factor=0.1)

    plotter = pv.Plotter(off_screen=args.out is not None)
    plotter.add_mesh(pv.Sphere(radius=1.0, theta_resolution=60, phi_resolution=60), color="lightgray", opacity=0.6)
    plotter.add_mesh(arrows, color="blue")
    plotter.show(screenshot=args.out, window_size=[1920, 1080])


if __name__ == "__main__":
    main()
