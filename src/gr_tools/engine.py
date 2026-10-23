"""Symbolic curvature engine: Christoffel symbols, Riemann, Ricci and Kretschmann."""
import sympy as sp
import itertools


class GREngine:
    """
    A symbolic General Relativity engine for computing spacetime geometry tensors.

    This class evaluates and caches the fundamental geometric quantities of a
    provided spacetime metric, including Christoffel symbols, the Riemann
    curvature tensor, the Ricci tensor, and the Kretschmann scalar.
    """

    def __init__(self, coords, metric):
        """
        Initializes the geometry engine.

        Args:
            coords (list): A list of SymPy symbols representing the coordinates.
            metric (sp.Matrix): The covariant metric tensor as a SymPy Matrix.
        """
        self.coords = coords
        self.g = metric
        self.n = len(coords)
        self.g_inv = metric.inv()

        self._Gamma = None
        self._Riemann = None
        self._Ricci = None
        self._Kretschmann = None

    def christoffel(self):
        """
        Computes the Christoffel symbols of the second kind.

        Returns:
            sp.MutableDenseNDimArray: A 3D array containing the Christoffel symbols.
        """
        if self._Gamma is not None:
            return self._Gamma

        n, coords, g, g_inv = self.n, self.coords, self.g, self.g_inv
        Gamma = sp.MutableDenseNDimArray.zeros(n, n, n)

        for rho, mu, nu in itertools.product(range(n), repeat=3):
            term = 0
            for sigma in range(n):
                if g_inv[rho, sigma] != 0:
                    term += g_inv[rho, sigma] * (
                            sp.diff(g[sigma, mu], coords[nu]) +
                            sp.diff(g[sigma, nu], coords[mu]) -
                            sp.diff(g[mu, nu], coords[sigma])
                    )
            if term != 0:
                Gamma[rho, mu, nu] = sp.simplify(term / 2)

        self._Gamma = Gamma
        return Gamma

    def riemann(self):
        """
        Computes the Riemann curvature tensor.

        Returns:
            sp.MutableDenseNDimArray: A 4D array containing the Riemann tensor components.
        """
        if self._Riemann is not None:
            return self._Riemann

        Gamma = self.christoffel()
        n, coords = self.n, self.coords
        Riem = sp.MutableDenseNDimArray.zeros(n, n, n, n)

        for rho, sigma, mu, nu in itertools.product(range(n), repeat=4):
            term = sp.diff(Gamma[rho, nu, sigma], coords[mu]) \
                   - sp.diff(Gamma[rho, mu, sigma], coords[nu])
            for lam in range(n):
                term += Gamma[rho, mu, lam] * Gamma[lam, nu, sigma]
                term -= Gamma[rho, nu, lam] * Gamma[lam, mu, sigma]
            if term != 0:
                Riem[rho, sigma, mu, nu] = sp.simplify(term)

        self._Riemann = Riem
        return Riem

    def ricci(self):
        """
        Computes the Ricci tensor by contracting the Riemann tensor.

        Returns:
            sp.MutableDenseMatrix: A 2D matrix containing the Ricci tensor components.
        """
        if self._Ricci is not None:
            return self._Ricci

        Riem = self.riemann()
        n = self.n
        Ric = sp.MutableDenseMatrix.zeros(n, n)

        for mu, nu in itertools.product(range(n), repeat=2):
            term = sum(Riem[rho, mu, rho, nu] for rho in range(n))
            if term != 0:
                Ric[mu, nu] = sp.simplify(term)

        self._Ricci = Ric
        return Ric

    def kretschmann(self):
        """
        Computes the Kretschmann scalar.

        The scalar is found by fully contracting the Riemann tensor with itself.
        Index raising and lowering is performed sequentially prior to contraction.

        Returns:
            sp.Expr: The simplified symbolic expression of the Kretschmann scalar.
        """
        if self._Kretschmann is not None:
            return self._Kretschmann

        R_uddd = self.riemann()
        n = self.n
        g, g_inv = self.g, self.g_inv

        R_dddd = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            term = sum(g[a, e] * R_uddd[e, b, c, d] for e in range(n) if g[a, e] != 0)
            if term != 0: R_dddd[a, b, c, d] = term

        R_uddd_raised = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            term = sum(g_inv[a, e] * R_dddd[e, b, c, d] for e in range(n) if g_inv[a, e] != 0)
            if term != 0: R_uddd_raised[a, b, c, d] = term

        R_uudd = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            term = sum(g_inv[b, e] * R_uddd_raised[a, e, c, d] for e in range(n) if g_inv[b, e] != 0)
            if term != 0: R_uudd[a, b, c, d] = term

        R_uuud = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            term = sum(g_inv[c, e] * R_uudd[a, b, e, d] for e in range(n) if g_inv[c, e] != 0)
            if term != 0: R_uuud[a, b, c, d] = term

        R_uuuu = sp.MutableDenseNDimArray.zeros(n, n, n, n)
        for a, b, c, d in itertools.product(range(n), repeat=4):
            term = sum(g_inv[d, e] * R_uuud[a, b, c, e] for e in range(n) if g_inv[d, e] != 0)
            if term != 0: R_uuuu[a, b, c, d] = term

        K = sum(R_dddd[a, b, c, d] * R_uuuu[a, b, c, d]
                for a, b, c, d in itertools.product(range(n), repeat=4))

        self._Kretschmann = sp.simplify(K)
        return self._Kretschmann
