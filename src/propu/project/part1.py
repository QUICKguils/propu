"""Code developed to answer the first part of the project.

Considering an engine rotational speed of 9000 RPM and an advance velocity of
20 m/s, compute the thrust, the power absorbed by the propeller and the
propulsive efficiency for each propeller.
Plot the thrust and power distributions along the span. Compare the propellers
and discuss in terms of collective pitch, velocity triangles, distributions of
torque and thrust.
"""

import matplotlib.pyplot as plt

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm


def main(*, sdiv=18, flag_out=True) -> dict[str, bemt.BemSolution]:
    """Execute the first part of the project."""
    sol_dict = dict()
    oper = get_operating_conditions()
    plot = sol_plotter(flag_out)
    display = sol_displayer(flag_out)

    for prop_name in stm.APC_dict:
        prop = stm.create_apc_propeller(prop_name)
        sol = bemt.bem(prop, oper, sdiv=sdiv)
        sol_dict[prop.keyword] = sol
        plot(sol)
        display(sol)

    return sol_dict


def get_operating_conditions() -> bemt.OperatingConditions:
    """Get the operating conditions, as specified in the statement."""
    Om = 9000 * uconv("rpm", "rad/s")  # Rotation speed [rad/s]
    v_inf = 20  # Wind speed [m/s]
    rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
    mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
    return bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)


def sol_displayer(flag_out=True):
    """Return a solution displayer, depending on `flag_out`."""

    def display(sol: bemt.BemSolution) -> None:
        """Display the solutions of the first part of the project."""
        print(
            f"Solution for {sol.prop.pretty_name}",
            f"    Thrust: {sol.thrust:.2f} N",
            f"    Power: {sol.power:.2f} W",
            f"    Efficiency: {(sol.eta*1e2):.2f} %",
            sep="\n"
        )

    def no_display(*args, **kwargs):
        pass

    if flag_out:
        return display
    return no_display


def sol_plotter(flag_out=True):
    """Return a solution plotter, depending on `flag_out`.

    The matplotlib object instance is kept as a closure,
    to ensure its idempotence through multiple plotter calls.
    """
    def plot(sol: bemt.BemSolution):
        """Plot the solutions of the first part of the project."""
        r_adim = sol.r_dist / sol.prop.geometry.span
        dTdr = sol.thrust_dist / sol.dr_dist
        dPdr = sol.power_dist / sol.dr_dist

        # Thrust distribution along the blades span
        thrust_line, = ax_thrust.plot(
            r_adim, dTdr,
            linewidth=0.8, label=sol.prop.pretty_name,
        )
        ax_thrust.scatter(
            r_adim[sol.flag_converged_dist], dTdr[sol.flag_converged_dist],
            s=30, marker=".", color=thrust_line.get_color(),
        )
        ax_thrust.scatter(
            r_adim[~sol.flag_converged_dist], dTdr[~sol.flag_converged_dist],
            s=30, marker="x", color=thrust_line.get_color(),
        )
        ax_thrust.set_title("Thrust distribution")
        ax_thrust.set_xlabel("$r/R$")
        ax_thrust.set_ylabel("$dT/dr$ [N/m]")
        ax_thrust.legend()

        # Power distribution along the blades span
        power_line, = ax_power.plot(
            r_adim, dPdr,
            linewidth=0.8, label=sol.prop.pretty_name,
        )
        ax_power.scatter(
            r_adim[sol.flag_converged_dist], dPdr[sol.flag_converged_dist],
            s=30, marker=".", color=power_line.get_color(),
        )
        ax_power.scatter(
            r_adim[~sol.flag_converged_dist], dPdr[~sol.flag_converged_dist],
            s=30, marker="x", color=power_line.get_color(),
        )
        ax_power.set_title("Power distribution")
        ax_power.set_xlabel("$r/R$")
        ax_power.set_ylabel("$dP/dr$ [W/m]")
        ax_power.legend()

        fig.show()

    def no_plot(*args, **kwargs):
        pass

    if flag_out:
        fig, (ax_thrust, ax_power) = plt.subplots(1, 2, figsize=(10, 5))
        return plot

    return no_plot
