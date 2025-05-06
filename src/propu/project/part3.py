"""Code developed to answer the third part of the project.

Consider a flight speed of 15 m/s and a drag of 11 N.
Determine the best propeller for this flight case and its rotational velocity.
What is the propulsive efficiency of the propeller in these conditions?
"""

from typing import NamedTuple

import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm


class Performance(NamedTuple):
    J: np.ndarray[float]
    eta: np.ndarray[float]
    drag: np.ndarray[float]
    converged: np.ndarray[bool]


class SolvedPoint(NamedTuple):
    J: float
    eta: float
    drag: float


class Solution(NamedTuple):
    """Associate a propeller to its computed performances."""
    prop: bemt.Propeller
    perf: Performance
    solved_point: SolvedPoint


def main(*, out_enabled=True) -> dict[str, Solution]:
    """Execute the third part of the project."""
    sol_dict = dict()

    for Propeller in stm.APCPropellers:
        prop = Propeller()
        perf = compute_performance(prop)
        solved_point = compute_solved_point(prop, perf)
        sol_dict[prop.keyword] = Solution(prop=prop, perf=perf, solved_point=solved_point)

    if out_enabled:
        plot_solution(sol_dict)

    return sol_dict


def get_operating_conditions(prop: bemt.Propeller, J: float) -> bemt.OperatingConditions:
    rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
    mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
    v_inf = 15  # Wind speed [m/s]
    D = 2 * prop.geometry.span
    Om = v_inf / J / D * uconv("rps", "rad/s")
    return bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)


def compute_performance(prop: bemt.Propeller) -> Performance:
    # NOTE:
    # APC recommends a max rpm of 150_000/D for the thin electric series,
    # where D is the diameter in inches.
    J_range = np.linspace(0.2, 1.2, 30)  # TODO: maybe pass this as argument
    eta = np.zeros(J_range.shape)
    drag = np.zeros(J_range.shape)
    converged = np.zeros(J_range.shape, dtype=bool)

    for i, J in enumerate(J_range):
        oper = get_operating_conditions(prop, J)
        bem_sol = bemt.bem(prop, oper)
        if bem_sol.eta < 0:  # Stop computations and trim array views
            J_range = J_range[:i]
            eta = eta[:i]
            drag = drag[:i]
            converged = converged[:i]
            break

        eta[i] = bem_sol.eta
        drag[i] = bem_sol.thrust
        converged[i] = np.all(bem_sol.converged_dist)

    return Performance(J=J_range, eta=eta, drag=drag, converged=converged)


def compute_solved_point(prop: bemt.Propeller, perf: Performance) -> SolvedPoint:
    DRAG_IMPOSED = 11  # [N]

    # First solution approximation
    # See this as a small minimization-like problem: f(x) = 0
    # -> solution correspond to the argmin of f(J) := drag(J) - DRAG_IMPOSED
    arg_aprox = np.argmin(np.abs(perf.drag - DRAG_IMPOSED))
    J = perf.J[arg_aprox]

    # More refined solution
    # Perform a simple bissection in the neighborhood of the first approx.
    J_bounds = [perf.J[arg_aprox-1], perf.J[arg_aprox+1]]
    max_iter = 10
    atol = 1e-4
    rtol = 1e-3
    has_converged = False
    for n_iter in range(1, max_iter + 1):
        oper = get_operating_conditions(prop, J)
        bem_sol = bemt.bem(prop, oper)
        drag = bem_sol.thrust

        # Stop the iterations if the drag is close enough to the imposed one
        if np.allclose(drag, DRAG_IMPOSED, atol=atol, rtol=rtol):
            has_converged = True
            break

        # Bissection update
        if drag > DRAG_IMPOSED:
            J_bounds[0] = J
        else:
            J_bounds[-1] = J
        J = np.mean(J_bounds)

    if not has_converged:
        print("yo seems the simple bissection is too wobbly wonky")

    return SolvedPoint(J=float(J), eta=float(bem_sol.eta), drag=float(drag))


def plot_solution(sol_dict: dict[str, Solution]) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    import matplotlib.pyplot as plt
    from propu.mplrc import REPORT_TW
    fig, (ax_eta, ax_drag) = plt.subplots(2, 1, figsize=(REPORT_TW, REPORT_TW))

    def plot(sol: Solution):
        """Plot the solutions of this project part."""
        (prop, perf, solved_point) = sol

        # Propulsive efficiencies
        eta_line, = ax_eta.plot(
            perf.J[perf.converged], perf.eta[perf.converged],
            linestyle="solid", linewidth=0.8,
            label=f"{prop.pretty_name}",
        )
        # ax_eta.scatter(
        #     perf.J[perf.converged], perf.eta[perf.converged],
        #     s=30, zorder=2.5, marker=".", color=eta_line.get_color(),
        # )
        # ax_eta.scatter(
        #     perf.J[~perf.converged], perf.eta[~perf.converged],
        #     s=30, zorder=2.5, marker="x", color=eta_line.get_color(),
        # )
        ax_eta.scatter(
            solved_point.J, solved_point.eta,
            s=30, zorder=2.5, marker="x", color=eta_line.get_color(),
        )
        ax_eta.set_xlabel(r"$J$")
        ax_eta.set_ylabel(r"$\eta_p$")

        # Total drag on the propeller
        drag_line, = ax_drag.plot(
            perf.J[perf.converged], perf.drag[perf.converged],
            linestyle="solid", linewidth=0.8,
        )
        # ax_drag.scatter(
        #     perf.J[perf.converged], perf.drag[perf.converged],
        #     s=30, zorder=2.5, marker=".", color=drag_line.get_color(),
        # )
        # ax_drag.scatter(
        #     perf.J[~perf.converged], perf.drag[~perf.converged],
        #     s=30, zorder=2.5, marker="x", color=drag_line.get_color(),
        # )
        ax_drag.scatter(
            solved_point.J, solved_point.drag,
            s=30, zorder=2.5, marker="x", color=drag_line.get_color(),
        )
        ax_drag.set_xlabel(r"$J$")
        ax_drag.set_ylabel(r"$D$ (N)")

    for _, sol in sol_dict.items():
        plot(sol)

    fig.legend(loc='outside upper center', ncols=4)
    fig.show()
