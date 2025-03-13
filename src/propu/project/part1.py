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

# Operating conditions
Om = 9000 * uconv("rpm", "rad/s")  # Rotation speed [rad/s]
v_inf = 20  # Wind speeds [m/s]
rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
oper = bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)


def main(*, sdiv=20, flag_plot=True):
    """Execute the first part of the project."""
    sol_dict = dict()
    plot = sol_plotter(flag_plot)

    for iprop, prop_name in enumerate(stm.APC_dict):
        prop = stm.create_apc_propeller(prop_name)
        sol = bemt.bem(prop, oper, sdiv=sdiv)
        sol_dict[prop.name] = sol
        plot(prop, sol)

    return sol_dict


def sol_plotter(flag_plot=True):
    """Return a solution plotter, depending on flag_plot.

    The matplotlib object instance is kept as closure,
    to ensure its idempotence.
    """
    fig, (ax_thrust, ax_power) = plt.subplots(1, 2, figsize=(10, 5))

    def plot(prop: bemt.Propeller, sol: bemt.BemSolution):
        """Plot the solutions of the first part of the project."""
        # Thrust distribution along the blades span
        ax_thrust.plot(
            sol.r_dist / prop.geometry.span,
            sol.thrust_dist / sol.dr_dist,
            marker=".", linewidth=0.5, label=prop.name,
        )
        ax_thrust.set_title("Thrust distribution")
        ax_thrust.set_xlabel("$r/R$")
        ax_thrust.set_ylabel("$dT/dr$ [N/m]")
        ax_thrust.legend()

        # Power distribution along the blades span
        ax_power.plot(
            sol.r_dist / prop.geometry.span,
            sol.torque_dist * Om / sol.dr_dist,
            marker=".", linewidth=0.5, label=prop.name,
        )
        ax_power.set_title("Power distribution")
        ax_power.set_xlabel("$r/R$")
        ax_power.set_ylabel("$dP/dr$ [W/m]")
        ax_power.legend()

        fig.show()

    def no_plot(*args, **kwargs):
        pass

    if flag_plot:
        return plot
    return no_plot
