"""Code developed to answer the second part of the project.

Plot the thrust and power coefficients, and propulsive efficiency against the
advance ratio. Compare your results with the experimental data from [1].
Discuss your results with the assumptions of the BEMT, flow conditions and
forces.
"""

from typing import NamedTuple
import re

import matplotlib.pyplot as plt
import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm

MEASUREMENTS = [
    # "apce_9x45_rd0996_4002.txt",
    # "apce_9x45_jb0998_6018.txt",
    # "apce_9x6_rd0988_4003.txt",
    # "apce_9x6_rd0991_6038.txt",
    "apce_11x7_kt0535_3003.txt",
    "apce_11x7_kt0539_4997.txt",
    # "apce_11x10_kt0511_3006.txt",
    # "apce_11x10_kt0463_5007.txt",
]


class Performance(NamedTuple):
    J: np.ndarray[float]
    CT: np.ndarray[float]
    CP: np.ndarray[float]
    eta: np.ndarray[float]


class Part2Solution(NamedTuple):
    """Gather all relevent input data and output results."""
    prop: bemt.Propeller
    Om: float
    measured: Performance
    computed: Performance
    converged: np.ndarray[bool]


def main(*, sdiv=18, flag_out=True) -> dict[str, Part2Solution]:
    """Execute the second part of the project."""
    sol_dict = dict()
    plot_sol = sol_plotter(flag_out)

    for fname in MEASUREMENTS:
        part2_sol = compute_performance(fname, sdiv=sdiv)
        sol_dict[fname] = part2_sol
        plot_sol(part2_sol)

    return sol_dict


def extract_apc_data(fname: str) -> tuple[str, float, Performance]:
    """Extract data from a UIUC performance measurement file."""
    # WARN: Good ol' regexp. Not robust, but works fine here.
    prop_name = re.search(r'([^_]+_[^_]+)', fname).group(1)
    rpm = float(re.search(r'_(\d+)\.txt$', fname).group(1))

    data = np.loadtxt(str(stm._DATA_PATH / "apc" / fname), skiprows=1)
    perf = Performance(J=data[:, 0], CT=data[:, 1], CP=data[:, 2], eta=data[:, 3])

    return (prop_name, rpm, perf)


def compute_performance(fname: str, *, sdiv) -> Part2Solution:
    (prop_name, rpm, perf_measured) = extract_apc_data(fname)
    prop = stm.create_apc_propeller(prop_name)

    rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
    mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
    n = rpm * uconv("rpm", "rps")
    Om = rpm * uconv("rpm", "rad/s")
    D = 2 * prop.geometry.span

    CT_computed = np.zeros(perf_measured.J.shape)
    CP_computed = np.zeros(perf_measured.J.shape)
    eta_computed = np.zeros(perf_measured.J.shape)
    converged_computed = np.zeros(perf_measured.J.shape, dtype=bool)

    for i, J in enumerate(perf_measured.J):
        v_inf = n * D * J
        oper = bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)
        bem_sol = bemt.bem(prop, oper, sdiv=sdiv)

        CT_computed[i] = bem_sol.thrust / (rho*n**2*D**4)
        CP_computed[i] = bem_sol.power / (rho*n**3*D**5)
        eta_computed[i] = J * CT_computed[i]/CP_computed[i]
        converged_computed[i] = np.all(bem_sol.flag_converged_dist)

    perf_computed = Performance(
        J=perf_measured.J,
        CT=CT_computed,
        CP=CP_computed,
        eta=eta_computed
    )

    return Part2Solution(
        prop=prop,
        Om=Om,
        measured=perf_measured,
        computed=perf_computed,
        converged=converged_computed
    )


# TODO: plot not converged sols with crosses
def sol_plotter(flag_out=True):
    """Return a solution plotter, depending on flag_plot.

    The matplotlib object instances are kept as a closures,
    to ensure its idempotence through multiple plotter calls.
    """
    def plot(sol: Part2Solution):
        """Plot the solutions of the first part of the project."""
        (prop, Om, measured, computed, converged) = sol
        rpm = Om * uconv("rad/s", "rpm")
        print(converged)

        # Thrust coefficients
        measured_CT_line, = axs["CT"].plot(
            measured.J, measured.CT,
            linestyle="dashed", linewidth=0.8,
        )
        axs["CT"].plot(
            computed.J, computed.CT,
            linestyle="solid", linewidth=0.8,
            color= measured_CT_line.get_color(),
        )
        axs["CT"].scatter(
            computed.J[converged], computed.CT[converged],
            s=30, marker=".", color=measured_CT_line.get_color(),
        )
        axs["CT"].scatter(
            computed.J[~converged], computed.CT[~converged],
            s=30, marker="x", color=measured_CT_line.get_color(),
        )

        axs["CT"].set_xlabel(r"$J$")
        axs["CT"].set_ylabel(r"$C_T$")

        # Power coefficients
        measured_CP_line, = axs["CP"].plot(
            measured.J, measured.CP,
            linestyle="dashed", linewidth=0.8,
        )
        axs["CP"].plot(
            computed.J, computed.CP,
            linestyle="solid", linewidth=0.8,
            color=measured_CP_line.get_color(),
        )
        axs["CP"].scatter(
            computed.J[converged], computed.CP[converged],
            s=30, marker=".", color=measured_CP_line.get_color(),
        )
        axs["CP"].scatter(
            computed.J[~converged], computed.CP[~converged],
            s=30, marker="x", color=measured_CP_line.get_color(),
        )
        axs["CP"].set_xlabel(r"$J$")
        axs["CP"].set_ylabel(r"$C_P$")

        # Propulsive efficiencies
        measured_eta_line, = axs["eta"].plot(
            measured.J, measured.eta,
            linestyle="dashed", linewidth=0.8,
        )
        axs["eta"].plot(
            computed.J, computed.eta,
            linestyle="solid", linewidth=0.8,
            color=measured_eta_line.get_color(),
            label=f"{prop.pretty_name} @ {rpm:.0f} rpm",
        )
        axs["eta"].scatter(
            computed.J[converged], computed.eta[converged],
            s=30, marker=".", color=measured_eta_line.get_color(),
        )
        axs["eta"].scatter(
            computed.J[~converged], computed.eta[~converged],
            s=30, marker="x", color=measured_eta_line.get_color(),
        )
        axs["eta"].set_xlabel(r"$J$")
        axs["eta"].set_ylabel(r"$\eta_p$")

        fig.legend(loc='outside upper center', ncols=2)
        fig.show()

    def no_plot(*args, **kwargs):
        pass

    if flag_out:
        from propu.mplrc import REPORT_TW
        fig, axs = plt.subplot_mosaic(
            [["CT", "CP"], ["eta", "eta"]],
            figsize=(REPORT_TW, REPORT_TW)
        )
        return plot

    return no_plot
