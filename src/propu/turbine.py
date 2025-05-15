"""Gas turbine analyses."""

from typing import NamedTuple
import warnings

import numpy as np

from propu import constant as c


def T_static2total(M, gamma=c.gamma_air):
    """Static to total temperature conversion factor."""
    return 1 + (gamma - 1) / 2 * M**2


def p_static2total(M, gamma=c.gamma_air):
    """Static to total pressure conversion factor."""
    return T_static2total(M, gamma) ** (gamma / (gamma - 1))


def mdot_chocking(pt, Tt, R, gamma, A):
    """Choking mass flow."""
    return pt / np.sqrt(gamma * R * Tt) * (2 / (gamma + 1)) ** ((gamma + 1) / 2 / (gamma - 1)) * A


def mdot_corrected(pt, Tt, m_flow):
    """Corrected/reduced mass flow."""
    return m_flow * (c.p_ref / pt) * np.sqrt(Tt / c.T_ref)


def bisection(f, bounds: tuple[float, float], *, max_iter=20, ftol=1e-3, xtol=1e-3):
    """Solve a minimization problem of the form `f(x) = 0` using the bisection method.

    This function finds a root of the function `f` within the specified bounds by iteratively
    narrowing down the interval where the root lies. The method assumes that `f` is monotonically
    increasing and possesses a root within the provided bounds.

    Parameters
    ----------
    f: function
        The function to minimize.
    bounds: tuple[float, float]
        The bounds for the searching range `x` values.
    max_iter: int, optional
        The maximum number of iterations that are tolerated.
    ftol : float, optional
        The absolute tolerance for the function value.
    xtol : float, optional
        The relative tolerance for the interval width.

    Returns
    -------
    Solution : NamedTuple
        A named tuple containing the following fields:
        - `x` : float
            The approximate root of the function.
        - `residual` : float
            The value of the function at the approximate root.
        - `converged` : bool
            Raised to True if `ftol` or `xtol` are satisfied.
        - `n_iter` : int
            The number of iterations performed.
    """
    class Solution(NamedTuple):
        x: float
        residual: float
        converged: bool
        n_iter: int

    converged = False
    for n_iter in range(1, max_iter + 1):
        x = (bounds[0] + bounds[-1]) / 2
        residual = f(x)

        # Check tolerances
        ftol_respected = ftol > abs(residual)
        xtol_respected = xtol > abs(bounds[0] - bounds[-1]) / max(abs(x), 1)
        if ftol_respected:
            break

        # Update bounds
        if residual > 0:
            bounds = (bounds[0], x)
        else:
            bounds = (x, bounds[-1])

    # Assess the convergence
    if ftol_respected:
        converged = True
    elif xtol_respected:
        converged = True
        warnings.warn("bisection -- xtol satisfied, but not ftol.")
    else:
        warnings.warn(f"bisection -- No convergence after {n_iter} iterations.")

    return Solution(x=x, residual=residual, converged=converged, n_iter=n_iter)
