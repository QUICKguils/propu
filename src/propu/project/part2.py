"""Code developed to answer the second part of the project.

Plot the thrust and power coefficients, and propulsive efficiency against the
advance ratio. Compare your results with the experimental data from [1].
Discuss your results with the assumptions of the BEMT, flow conditions and
forces.
"""

import re
from typing import NamedTuple

import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.project import statement as stm

MEASUREMENTS = [
    "apce_9x45_rd0996_4002.txt",
    # "apce_9x45_jb0998_6018.txt",
    "apce_9x6_rd0988_4003.txt",
    # "apce_9x6_rd0991_6038.txt",
    "apce_11x7_kt0535_3003.txt",
    # "apce_11x7_kt0539_4997.txt",
    "apce_11x10_kt0511_3006.txt",
    # "apce_11x10_kt0463_5007.txt",
]


class Performance(NamedTuple):
    J: np.ndarray[float]
    CT: np.ndarray[float]
    CP: np.ndarray[float]
    eta: np.ndarray[float]


class Solution(NamedTuple):
    prop: bemt.Propeller
    Om: float
    measured: Performance
    computed: Performance
    converged: np.ndarray[bool]


def main(*, out_enabled=True) -> dict[str, Solution]:
    """Execute the second part of the project."""
    sol_dict = dict()

    for fname in MEASUREMENTS:
        sol_dict[fname] = compute_performance(fname)

    if out_enabled:
        plot_solution(sol_dict)

    return sol_dict


def extract_apc_data(fname: str) -> tuple[str, float, Performance]:
    """Extract data from a UIUC performance measurement file."""
    # WARN: Good ol' regexp. Not robust, but works fine here.
    prop_name = re.search(r"([^_]+_[^_]+)", fname).group(1)
    rpm = float(re.search(r"_(\d+)\.txt$", fname).group(1))

    data = np.loadtxt(str(stm._DATA_PATH / "apc" / fname), skiprows=1)
    perf = Performance(J=data[:, 0], CT=data[:, 1], CP=data[:, 2], eta=data[:, 3])

    return (prop_name, rpm, perf)


def compute_performance(fname: str) -> Solution:
    (prop_keyword, rpm, perf_measured) = extract_apc_data(fname)
    prop = stm.create_apc_propeller(prop_keyword)

    n = rpm * uconv("rpm", "rps")
    Om = rpm * uconv("rpm", "rad/s")
    D = 2 * prop.geometry.span

    CT = np.zeros(perf_measured.J.shape)
    CP = np.zeros(perf_measured.J.shape)
    eta = np.zeros(perf_measured.J.shape)
    converged = np.zeros(perf_measured.J.shape, dtype=bool)

    for i, J in enumerate(perf_measured.J):
        v_inf = n * D * J
        oper = bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=stm.rho, mu=stm.mu)
        bem_sol = bemt.bem(prop, oper)

        CT[i] = bem_sol.thrust / (stm.rho * n**2 * D**4)
        CP[i] = bem_sol.power / (stm.rho * n**3 * D**5)
        eta[i] = bem_sol.eta  # alternatively, eta[i] = J * CT[i]/CP[i]
        converged[i] = np.all(bem_sol.converged_dist)

    perf_computed = Performance(J=perf_measured.J, CT=CT, CP=CP, eta=eta)

    return Solution(
        prop=prop, Om=Om, measured=perf_measured, computed=perf_computed, converged=converged
    )


def plot_solution(sol_dict: dict[str, Solution]) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    import matplotlib.pyplot as plt

    from propu.mplrc import REPORT_TW

    fig, axs = plt.subplot_mosaic([["CT", "CP"], ["eta", "eta"]], figsize=(REPORT_TW, REPORT_TW))

    def plot(sol: Solution):
        """Plot the solutions of this project part."""
        (prop, Om, measured, computed, converged) = sol
        rpm = Om * uconv("rad/s", "rpm")

        # Thrust coefficients
        (computed_CT_line,) = axs["CT"].plot(
            computed.J, computed.CT,
            linestyle="solid", linewidth=0.8,
        )
        axs["CT"].scatter(
            computed.J[converged], computed.CT[converged],
            s=15, linewidths=1, zorder=2.5, marker=".", color=computed_CT_line.get_color(),
        )
        axs["CT"].scatter(
            computed.J[~converged], computed.CT[~converged],
            s=15, linewidths=1, zorder=2.5, marker="x", color=computed_CT_line.get_color(),
        )
        axs["CT"].scatter(
            measured.J, measured.CT,
            s=15, linewidths=1, zorder=2.5, marker="+", color=computed_CT_line.get_color(),
        )

        axs["CT"].set_xlabel(r"$J$")
        axs["CT"].set_ylabel(r"$C_T$")

        # Power coefficients
        (computed_CP_line,) = axs["CP"].plot(
            computed.J, computed.CP,
            linestyle="solid", linewidth=0.8,
        )
        axs["CP"].scatter(
            computed.J[converged], computed.CP[converged],
            s=15, linewidths=1, zorder=2.5, marker=".", color=computed_CP_line.get_color(),
        )
        axs["CP"].scatter(
            computed.J[~converged], computed.CP[~converged],
            s=15, linewidths=1, zorder=2.5, marker="x", color=computed_CP_line.get_color(),
        )
        axs["CP"].scatter(
            measured.J, measured.CP,
            s=15, linewidths=1, zorder=2.5, marker="+", color=computed_CP_line.get_color(),
        )
        axs["CP"].set_xlabel(r"$J$")
        axs["CP"].set_ylabel(r"$C_P$")

        # Propulsive efficiencies
        (computed_eta_line,) = axs["eta"].plot(
            computed.J, computed.eta,
            linestyle="solid", linewidth=0.8, label=f"{prop.pretty_name} @ {rpm:.0f} rpm",
        )
        axs["eta"].scatter(
            computed.J[converged], computed.eta[converged],
            s=15, linewidths=1, zorder=2.5, marker=".", color=computed_eta_line.get_color(),
        )
        axs["eta"].scatter(
            computed.J[~converged], computed.eta[~converged],
            s=15, linewidths=1, zorder=2.5, marker="x", color=computed_eta_line.get_color(),
        )
        axs["eta"].scatter(
            measured.J[measured.eta > 0], measured.eta[measured.eta > 0],
            s=15, linewidths=1, zorder=2.5, marker="+", color=computed_eta_line.get_color(),
        )
        axs["eta"].set_xlabel(r"$J$")
        axs["eta"].set_ylabel(r"$\eta_p$")

    for _, sol in sol_dict.items():
        plot(sol)

    fig.legend(loc="outside upper center", ncols=2)
    fig.show()
