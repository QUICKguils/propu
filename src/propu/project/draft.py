"""Code developed to answer the fourth part of the project.

For the latter configuration, investigate strategies to maintain a constant
rotational velocity of 8000 RPM with a fixed power input of 0.31 hp and varying
advance velocity.
"""

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray

from propu import bemt
from propu.constant import uconv
from propu.project import statement as stm

OM_IMPOSED = 8000 * uconv("rpm", "rad/s")
P_IMPOSED = 0.31 * uconv("hp", "W")


class Performance(NamedTuple):
    J: NDArray[np.float64]
    power: NDArray[np.float64]
    converged: NDArray[np.bool]


class SolvedPoint(NamedTuple):
    J: float
    power: float
    v_inf: float


class Solution(NamedTuple):
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
        display_solution(sol_dict)
        plot_solution(sol_dict)

    return sol_dict


def get_operating_conditions(prop: bemt.Propeller, J: float) -> bemt.OperatingConditions:
    D = 2 * prop.geometry.span
    v_inf = J * D * OM_IMPOSED * uconv("rad/s", "rps")
    return bemt.OperatingConditions(Om=OM_IMPOSED, v_inf=v_inf, rho=stm.rho, mu=stm.mu)


def compute_performance(prop: bemt.Propeller, J_range=np.linspace(0.2, 1, 30)) -> Performance:
    # NOTE:
    # APC recommends a max rpm of 150_000/D for the thin electric series,
    # where D is the diameter in inches.
    power = np.zeros(J_range.shape)
    converged = np.zeros(J_range.shape, dtype=bool)

    for i, J in enumerate(J_range):
        oper = get_operating_conditions(prop, J)
        bem_sol = bemt.bem(prop, oper)
        if bem_sol.eta < 0:  # Stop computations and trim array views
            J_range = J_range[:i]
            power = power[:i]
            converged = converged[:i]
            break

        power[i] = bem_sol.power
        converged[i] = np.all(bem_sol.converged_dist)

    return Performance(J=J_range, power=power, converged=converged)


def compute_solved_point(prop: bemt.Propeller, perf: Performance) :
    # First solution approximation
    # See this as a small minimization-like problem: f(x) = 0
    # -> solution correspond to the argmin of f(J) := power(J) - P_IMPOSED
    arg_aprox = np.argmin(np.abs(perf.power - P_IMPOSED))
    J = perf.J[arg_aprox]

    # More refined solution
    # Perform a simple bissection in the neighborhood of the first approx.
    J_bounds = [perf.J[arg_aprox - 1], perf.J[arg_aprox + 1]]  # FIX: can be out of array bounds
    max_iter = 10
    atol = 1e-2
    rtol = 1e-3
    has_converged = False
    for n_iter in range(1, max_iter + 1):
        oper = get_operating_conditions(prop, J)
        bem_sol = bemt.bem(prop, oper)
        power = bem_sol.power

        # Stop the iterations if the drag is close enough to the imposed one
        if np.allclose(power, P_IMPOSED, atol=atol, rtol=rtol):
            has_converged = True
            break

        # Bissection update
        if power > P_IMPOSED:
            J_bounds[0] = J
        else:
            J_bounds[-1] = J
        J = np.mean(J_bounds)

    if not has_converged:
        print("yo seems the simple bissection is too wobbly wonky")

    D = 2 * prop.geometry.span
    v_inf = J * D * OM_IMPOSED * uconv("rad/s", "rps")

    return SolvedPoint(J=float(J), power=float(power), v_inf=v_inf)


def display_solution(sol_dict: dict[str, Solution]) -> None:
    """Pretty-print the computed solutions."""
    print("Part4 solution")
    for _, sol in sol_dict.items():
        print(
            f"Propeller: {sol.prop.pretty_name}",
            f"    Imposed rot. speed: {OM_IMPOSED * uconv("rad/s", "rpm"):.0f} rpm",
            f"    Imposed power: {P_IMPOSED:.2f} W",
            f"    Advance ratio: {sol.solved_point.J}",
            f"    Upstream wind speed: {sol.solved_point.v_inf} m/s",
            sep="\n",
        )


def plot_solution(sol_dict: dict[str, Solution]) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    import matplotlib.pyplot as plt
    from propu.mplrc import REPORT_TW

    fig, ax = plt.subplots()

    def plot(sol: Solution):
        (prop, perf, solved_point) = sol

        (line,) = ax.plot(
            perf.J[perf.converged], perf.power[perf.converged],
            linestyle="solid", linewidth=0.8, label=f"{prop.pretty_name}",
        )
        ax.scatter(
            perf.J[~perf.converged], perf.power[~perf.converged],
            s=20, linewidths=1, zorder=2.5, marker="x", color=line.get_color(),
        )
        ax.scatter(
            solved_point.J, solved_point.power,
            s=20, linewidths=1, zorder=2.5, marker="+", color=line.get_color(),
        )
        ax.set_xlabel(r"$J$")
        ax.set_ylabel(r"$P$ (W)")

    for _, sol in sol_dict.items():
        plot(sol)

    fig.legend(loc="outside upper center", ncols=4)
    fig.show()
