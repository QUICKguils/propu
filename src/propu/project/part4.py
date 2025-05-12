"""Code developed to answer the fourth part of the project.

For the latter configuration, investigate strategies to maintain a constant rotational velocity of
8000 RPM with a fixed power input of 0.31 hp and varying advance velocity.
"""

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from propu import bemt
from propu.constant import uconv
from propu.project import statement as stm

OM_IMPOSED = 8000 * uconv("rpm", "rad/s")
P_IMPOSED = 0.31 * uconv("hp", "W")
KEYWORD_PROP_IMPOSED = "apce_11x7"  # As determined in part3


def create_variable_APC(theta: float) -> bemt.Propeller:
    """Create an APC propeller, with fictious collective pitch.

    Parameters
    ----------
    theta: float
        The collective pitch angle, in radians.

    Notes
    -----
    The function uses the APC propeller already instantiated in the statement,
    to avoid the heavy computations of reading each time the data files in stm._DATA_PATH.
    """
    apc_11x7 = stm.get_apc_propeller(KEYWORD_PROP_IMPOSED)
    geometry = bemt.PropellerGeometry(
        span=apc_11x7.geometry.span,
        n_blades=apc_11x7.geometry.n_blades,
        stations=apc_11x7.geometry.stations,
        chords=apc_11x7.geometry.chords,
        pitches=apc_11x7.geometry.pitches + theta,
    )
    return bemt.Propeller(
        keyword=apc_11x7.keyword,
        pretty_name=apc_11x7.pretty_name,
        airfoil=apc_11x7.airfoil,
        geometry=geometry,
    )


class Point(NamedTuple):
    J: float
    power: float
    theta: float
    converged: bool


class Curve(NamedTuple):
    J: NDArray[np.float64]
    power: NDArray[np.float64]
    theta: NDArray[np.float64]
    converged: NDArray[np.bool]


class Contour(NamedTuple):
    J: NDArray[np.float64]
    theta: NDArray[np.float64]
    power: NDArray[np.float64]


class Solution(NamedTuple):
    curve: Curve
    contour: Contour


def main(*, out_enabled=True) -> Solution:
    """Execute the fourth part of the project."""
    J_range = np.linspace(0.2, 1.0, 20)
    thetad_range = np.linspace(0, 10, 15)
    theta_range = np.deg2rad(thetad_range)

    curve_sol = compute_curve(J_range)
    contour_sol = compute_contour(J_range, theta_range)
    sol = Solution(curve=curve_sol, contour=contour_sol)

    if out_enabled:
        plot_solution(sol)

    return sol


def get_operating_conditions(J: float) -> bemt.OperatingConditions:
    D = 2 * stm.get_apc_propeller(KEYWORD_PROP_IMPOSED).geometry.span
    v_inf = J * D * OM_IMPOSED * uconv("rad/s", "rps")
    return bemt.OperatingConditions(Om=OM_IMPOSED, v_inf=v_inf, rho=stm.rho, mu=stm.mu)


def compute_point(J: float, thetad_bounds=[0, 20]) -> Point:
    theta_bounds = np.deg2rad(thetad_bounds)

    # Perform a secant search in the neighborhood of the initial guess.
    # See this as a small minimization-like problem:
    #   f(theta) = 0, where f(theta) := power(theta) - P_IMPOSED.
    max_iter = 15
    rtol = 1e-3
    converged = False
    for _ in range(1, max_iter + 1):
        theta = np.mean(theta_bounds)
        prop = create_variable_APC(theta)
        oper = get_operating_conditions(J)
        power = bemt.bem(prop, oper).power

        # Stop the iterations if the power is close enough to the imposed one
        if np.allclose(power, P_IMPOSED, rtol=rtol):
            converged = True
            break

        # Bissection update
        if power > P_IMPOSED:
            theta_bounds[-1] = theta
        else:
            theta_bounds[0] = theta

    return Point(J=J, power=power, theta=theta, converged=converged)


def compute_curve(J_range: NDArray) -> Curve:
    power_range = np.zeros(J_range.shape)
    theta_range = np.zeros(J_range.shape)
    converged = np.zeros(J_range.shape, dtype=np.bool)

    for i, J in enumerate(J_range):
        solved_point = compute_point(J)
        power_range[i] = solved_point.power
        theta_range[i] = solved_point.theta
        converged[i] = solved_point.converged

    return Curve(J=J_range, power=power_range, theta=theta_range, converged=converged)


def compute_contour(J_range: NDArray, theta_range: NDArray) -> Contour:
    (J_grid, theta_grid) = np.meshgrid(J_range, theta_range)
    power_grid = np.zeros(J_grid.shape, dtype=np.float64)

    for i, (J, theta) in enumerate(zip(J_grid.flat, theta_grid.flat)):
        oper = get_operating_conditions(J)
        prop = create_variable_APC(theta)
        power_grid.flat[i] = bemt.bem(prop, oper).power

    return Contour(J=J_grid, theta=theta_grid, power=power_grid)


def plot_solution(sol: Solution) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    import matplotlib.pyplot as plt

    from propu.mplrc import REPORT_TW

    fig, ax = plt.subplots(figsize=(0.5 * REPORT_TW, 0.5 * REPORT_TW))

    topo = ax.contour(sol.contour.J, np.rad2deg(sol.contour.theta), sol.contour.power)
    (line,) = ax.plot(sol.curve.J, np.rad2deg(sol.curve.theta), linestyle="solid", linewidth=0.8)
    ax.scatter(
        sol.curve.J[sol.curve.converged], np.rad2deg(sol.curve.theta[sol.curve.converged]),
        s=15, linewidths=1, zorder=2.5, marker=".", color=line.get_color(),
    )
    ax.scatter(
        sol.curve.J[~sol.curve.converged], np.rad2deg(sol.curve.theta[~sol.curve.converged]),
        s=15, linewidths=1, zorder=2.5, marker="x", color=line.get_color(),
    )

    ax.clabel(topo, fontsize=9, fmt="%3.0f W")
    ax.set_xlabel(r"$J$")
    ax.set_ylabel(r"$\theta/\degree$")

    fig.show()


def _get_v_inf(sol: Solution) -> NDArray:
    D = 2 * stm.get_apc_propeller(KEYWORD_PROP_IMPOSED).geometry.span
    return sol.curve.J * D * OM_IMPOSED * uconv("rad/s", "rps")
