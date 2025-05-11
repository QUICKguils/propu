"""Answer the first part of the project.

Considering an engine rotational speed of 9000 RPM and an advance velocity of 20 m/s, compute the
thrust, the power absorbed by the propeller and the propulsive efficiency for each propeller.
Plot the thrust and power distributions along the span.
Compare the propellers and discuss in terms of collective pitch, velocity triangles, distributions
of torque and thrust.
"""

import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.project import statement as stm

V_IMPOSED = 20  # [m/s]
OM_IMPOSED = 9000 * uconv("rpm", "rad/s")


def main(*, out_enabled=True) -> dict[str, bemt.BemSolution]:
    """Execute the first part of the project."""
    sol_dict: dict[str, bemt.BemSolution] = dict()
    oper = bemt.OperatingConditions(Om=OM_IMPOSED, v_inf=V_IMPOSED, rho=stm.rho, mu=stm.mu)

    for prop in stm.propellers:
        sol_dict[prop.keyword] = bemt.bem(prop, oper)

    if out_enabled:
        display_solution(sol_dict)
        plot_solution(sol_dict)

    return sol_dict


def display_solution(sol_dict: dict[str, bemt.BemSolution]) -> None:
    """Pretty-print the computed solutions."""
    print(
        "Part 1 solution"
        f" (Rot. speed: {OM_IMPOSED * uconv('rad/s', 'rpm'):.0f} rpm,"
        f" v_inf: {V_IMPOSED} m/s)"
    )
    for _, sol in sol_dict.items():
        print(
            f"Propeller: {sol.prop.pretty_name}",
            f"    Thrust: {sol.thrust:.2f} N",
            f"    Power: {sol.power:.2f} W",
            f"    Efficiency: {(sol.eta * 1e2):.2f} %",
            sep="\n",
        )


def plot_solution(sol_dict: dict[str, bemt.BemSolution]) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    import matplotlib.pyplot as plt

    fig, (ax_thrust, ax_power) = plt.subplots(1, 2)

    def plot(sol: bemt.BemSolution) -> None:
        """Plot one computed solution."""
        r_adim = sol.r_dist / sol.prop.geometry.span
        dTdr = sol.thrust_dist / sol.dr_dist
        dPdr = sol.power_dist / sol.dr_dist

        # Thrust distribution along the blades span
        (thrust_line,) = ax_thrust.plot(
            r_adim, dTdr,
            linewidth=0.8, label=sol.prop.pretty_name,
        )
        ax_thrust.scatter(
            r_adim[sol.converged_dist], dTdr[sol.converged_dist],
            s=15, linewidths=1, zorder=2.5, marker=".", color=thrust_line.get_color(),
        )
        ax_thrust.scatter(
            r_adim[~sol.converged_dist], dTdr[~sol.converged_dist],
            s=15, linewidths=1, zorder=2.5, marker="x", color=thrust_line.get_color(),
        )
        ax_thrust.set_xlabel("$r/R$")
        ax_thrust.set_ylabel("$dT/dr$ [N/m]")

        # Power distribution along the blades span
        (power_line,) = ax_power.plot(
            r_adim, dPdr,
            linewidth=0.8,
        )
        ax_power.scatter(
            r_adim[sol.converged_dist], dPdr[sol.converged_dist],
            s=15, linewidths=1, zorder=2.5, marker=".", color=power_line.get_color(),
        )
        ax_power.scatter(
            r_adim[~sol.converged_dist], dPdr[~sol.converged_dist],
            s=15, linewidths=1, zorder=2.5, marker="x", color=power_line.get_color(),
        )
        ax_power.set_xlabel("$r/R$")
        ax_power.set_ylabel("$dP/dr$ [W/m]")

    for _, sol in sol_dict.items():
        plot(sol)

    fig.legend(loc="outside upper center", ncols=4)
    fig.show()


def plot_pitches() -> None:
    import matplotlib.pyplot as plt
    from propu.mplrc import REPORT_TW

    def find_optimal_aoa() -> float:
        arg_aoa = int(np.mean(np.argmax(stm.airfoil.cl / stm.airfoil.cd, axis=0)))
        # arg_aoa = np.argmin(stm.airfoil.cd, axis=0)[0]
        return stm.airfoil.aoa[arg_aoa]

    def plot(prop, ax):
        ax.plot(
            prop.geometry.stations/prop.geometry.span,
            np.rad2deg(prop.geometry.pitches),
            linewidth=0.8, label="Blade pitch"
        )
        ax.plot(
            sol.r_dist/prop.geometry.span,
            np.rad2deg(np.pi/2 + aoa_opt + sol.beta_dist),
            linewidth=0.8, label="Optimal pitch"
        )
        ax.set_title(f"{prop.pretty_name}")
        ax.set_xlabel(r"$r/R$")
        ax.set_ylabel(r"$\bar{\chi}$")

    oper = bemt.OperatingConditions(Om=OM_IMPOSED, v_inf=V_IMPOSED, rho=stm.rho, mu=stm.mu)
    aoa_opt = find_optimal_aoa()

    fig, axs = plt.subplots(2, 2, figsize=(REPORT_TW, REPORT_TW))

    for i, prop in enumerate(stm.propellers):
        sol = bemt.bem(prop, oper)
        plot(prop, axs.flat[i])

    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncols=2)
    fig.show()


