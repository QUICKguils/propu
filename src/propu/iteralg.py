"""Iteration algorithm utilities."""

from typing import NamedTuple
from warnings import warn


class IterTable:
    """Iteration table printer."""

    def __init__(self, variables: tuple, units: tuple):
        self.variables = variables
        self.units = units

        self.colspaces = tuple(
            max(10, len(v), len(u)) + 2 for (v, u) in zip(self.variables, self.units)
        )
        self.hline = "-" * (sum(self.colspaces) + len(self.variables) + 1)
        self.variables_row = "|" + "".join(
            f"{var:^{cs}}|" for (var, cs) in zip(self.variables, self.colspaces)
        )
        self.units_row = "|" + "".join(
            f"{unit:^{cs}}|" for (unit, cs) in zip(self.units, self.colspaces)
        )
        self.header = (self.hline, self.variables_row, self.units_row, self.hline)
        self.footer = (self.hline,)

        self.rows = []  # The future rows to add

    def add_row(self, *variables_value):
        if len(variables_value) is not len(self.variables):
            warn("Bro use my printing table correctly pls")
        self.rows.append(
            "|" + "".join(f"{val:^{cs}.4g}|" for (val, cs) in zip(variables_value, self.colspaces))
        )

    def print(self):
        print(*self.header, sep="\n")
        print(*self.rows, sep="\n")
        print(*self.footer)


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

    Raises
    ------
    ValueError
        If f(x) at bounds are not of opposite signs.
    """

    class Solution(NamedTuple):
        x: float
        residual: float
        converged: bool
        n_iter: int

    if f(bounds[0]) * f(bounds[-1]) > 0:
        raise ValueError("Seems the root is not contained in the given `bounds`")

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
