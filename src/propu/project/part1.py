"""Answer the first part of the project.

Considering an engine rotational speed of 9000 RPM and an advance velocity of
20 m/s, compute the thrust, the power absorbed by the propeller and the
propulsive efficiency for each propeller.
Plot the thrust and power distributions along the span. Compare the propellers
and discuss in terms of collective pitch, velocity triangles, distributions of
torque and thrust.
"""

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm


def main(*, out_enabled=True) -> dict[str, bemt.BemSolution]:
    """Execute the first part of the project."""
    sol_dict = dict()
    oper = get_operating_conditions()

    for Propeller in stm.APCPropellers:
        prop = Propeller()
        sol_dict[prop.keyword] = bemt.bem(prop, oper)

    if out_enabled:
        display_solution(sol_dict)
        plot_solution(sol_dict)

    return sol_dict


def get_operating_conditions() -> bemt.OperatingConditions:
    """Get the operating conditions, as specified in the statement."""
    Om = 9000 * uconv("rpm", "rad/s")  # Rotation speed [rad/s]
    v_inf = 20  # Wind speed [m/s]
    rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
    mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
    return bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)


def display_solution(sol_dict: dict[str, bemt.BemSolution]) -> None:
    """Pretty-print the computed solutions."""
    for _, sol in sol_dict.items():
        print(
            f"Solution for {sol.prop.pretty_name}",
            f"    Thrust: {sol.thrust:.2f} N",
            f"    Power: {sol.power:.2f} W",
            f"    Efficiency: {(sol.eta*1e2):.2f} %",
            sep="\n"
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
        thrust_line, = ax_thrust.plot(
            r_adim, dTdr,
            linewidth=0.8, label=sol.prop.pretty_name,
        )
        ax_thrust.scatter(
            r_adim[sol.converged_dist], dTdr[sol.converged_dist],
            s=30, zorder=2.5, marker=".", color=thrust_line.get_color(),
        )
        ax_thrust.scatter(
            r_adim[~sol.converged_dist], dTdr[~sol.converged_dist],
            s=30, zorder=2.5, marker="x", color=thrust_line.get_color(),
        )
        ax_thrust.set_xlabel("$r/R$")
        ax_thrust.set_ylabel("$dT/dr$ [N/m]")

        # Power distribution along the blades span
        power_line, = ax_power.plot(
            r_adim, dPdr,
            linewidth=0.8,
        )
        ax_power.scatter(
            r_adim[sol.converged_dist], dPdr[sol.converged_dist],
            s=30, zorder=2.5, marker=".", color=power_line.get_color(),
        )
        ax_power.scatter(
            r_adim[~sol.converged_dist], dPdr[~sol.converged_dist],
            s=30, zorder=2.5, marker="x", color=power_line.get_color(),
        )
        ax_power.set_xlabel("$r/R$")
        ax_power.set_ylabel("$dP/dr$ [W/m]")

    for _, sol in sol_dict.items():
        plot(sol)

    fig.legend(loc='outside upper center', ncols=4)
    fig.show()
