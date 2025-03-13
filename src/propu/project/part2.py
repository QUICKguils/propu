"""Code developed to answer the second part of the project.

Plot the thrust and power coefficients, and propulsive efficiency against the
advance ratio. Compare your results with the experimental data from [1].
Discuss your results with the assumptions of the BEMT, flow conditions and
forces.
"""

import matplotlib.pyplot as plt
import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm

# Operating conditions
Om_range = np.array([1000, 1000, 800, 800, 800, 700, 700]) * uconv("rpm", "rad/s")  # Rotation speeds [rad/s]
theta_75deg_range = np.array([15, 20, 25, 30, 35, 40, 45])  # Collective pitches [°]
theta_75_range = theta_75deg_range * uconv("deg", "rad")  # Collective pitches [rad]
v_inf_range = np.linspace(35, 115, 20) * uconv("mi/hr", "m/s")  # Wind speeds [m/s]
rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]


def main(*, sdiv=20, plot=True):
    """Execute the second part of the project."""
    thrusts, powers, Oms, v_infs = reproduce_report(sdiv=sdiv)
    J, C_T, C_P, eta = performance_coefficients(thrusts, powers, Oms, v_infs)

    # Plot the quantities, if desired
    plot and plot_sol(J, C_T, C_P, eta)

    return J, C_T, C_P, eta


# NOTE:
# Inspired from the simple indirect search,
# as used in the Aeroelasticity course.
def reproduce_report(*, sdiv=20):
    """Reproduce the experimental method conducted in the reference NACA report.

    Calculate the propeller thrust and power for an increasing wind speed, up
    to a maximum of 115 mph. Then decrease the propeller rotation speed to
    further increase the advance ratio. Stop the computation when the computed
    thrust becomes negative.
    """
    # Need for python list to be able to append dynamically, but
    # need for numpy array for array computation.
    thrusts = [np.array([]) for _ in theta_75_range]  # [N]
    powers = [np.array([]) for _ in theta_75_range]  # [W]
    Oms = [np.array([]) for _ in theta_75_range]  # [rad/s]
    v_infs = [np.array([]) for _ in theta_75_range]  # [m/s]

    for theta_id, _ in enumerate(theta_75_range):
        Om = Om_range[theta_id]

        # Loop over all prescribed wind speeds.
        for v_inf_id, v_inf in enumerate(v_inf_range):
            oper = bemt.OperatingConditions(
                Om=Om, theta_75=theta_75_range[theta_id], v_inf=v_inf, rho=rho, mu=mu
            )
            sol = bemt.bem(stm.prop, oper, sdiv=sdiv)

            if sol.thrust < 0:
                break

            Oms[theta_id] = np.append(Oms[theta_id], oper.Om)
            v_infs[theta_id] = np.append(v_infs[theta_id], oper.v_inf)
            thrusts[theta_id] = np.append(thrusts[theta_id], sol.thrust)
            powers[theta_id] = np.append(powers[theta_id], sol.torque * Om)

        # Algorithm heuristics
        dOm = 0.03 * Om
        dOm_refine_factor = 2
        thrust_atol = 1e-1
        max_iter = 40

        # Decrease the rotation speed to reach zero thrust.
        has_converged = False
        for n_iter in range(1, max_iter + 1):
            oper = bemt.OperatingConditions(
                Om=Om, theta_75=theta_75_range[theta_id], v_inf=v_inf, rho=rho, mu=mu
            )
            sol = bemt.bem(stm.prop, oper, sdiv=sdiv)

            if sol.thrust > 0:
                Oms[theta_id] = np.append(Oms[theta_id], oper.Om)
                v_infs[theta_id] = np.append(v_infs[theta_id], oper.v_inf)
                thrusts[theta_id] = np.append(thrusts[theta_id], sol.thrust)
                powers[theta_id] = np.append(powers[theta_id], sol.torque * Om)

            # Overshoot the thrust zero value ?
            # -> Step back and refine the Om step size.
            if sol.thrust < 0:
                Om += dOm
                dOm /= dOm_refine_factor
            else:
                Om -= dOm

            # Obtain convergence ?
            if np.abs(sol.thrust) < thrust_atol:
                has_converged = True
                break

        if not has_converged:
            print(f"project.part2 -- No convergence after {n_iter} iterations.")

    return thrusts, powers, Oms, v_infs


def performance_coefficients(thrusts, powers, Oms, v_infs):
    """Compute the performance coefficients."""

    J = [np.array([]) for _ in theta_75_range]  # Advance ratios
    C_T = [np.array([]) for _ in theta_75_range]  # Thrust coefficients
    C_P = [np.array([]) for _ in theta_75_range]  # Power coefficients
    eta = [np.array([]) for _ in theta_75_range]  # Propulsive efficiencies

    for i, _ in enumerate(theta_75_range):
        n = Oms[i] * uconv("rad/s", "rps")
        J[i] = v_infs[i] / (2 * stm.span * n)
        C_T[i] = 4 * thrusts[i] / ((2 * stm.span) ** 4 * rho * n**2)
        C_P[i] = 4 * powers[i] / ((2 * stm.span) ** 5 * rho * n**3)
        eta[i] = C_T[i] / C_P[i] * J[i]

    # NOTE: to translate in imperial units, use:
    # J   = J   / C.mph2ms * C.ft2m * C.tr2rad
    # C_T = C_T / C.lbf2n * C.ft2m**4 * C.slug2kg / C.ft2m**3 * C.tr2rad**2
    # C_P = C_P / C.hp2w * C.ft2m**5 * C.slug2kg / C.ft2m**3 * C.tr2rad**3

    return J, C_T, C_P, eta


def plot_sol(J, C_T, C_P, eta) -> None:
    """Plot the solutions of the second part of the project."""
    from propu.mplrc import REPORT_TW

    fig, axs = plt.subplot_mosaic(
        [["eta", "eta"], ["thrust", "power"]],
        figsize=(REPORT_TW, 0.8 * REPORT_TW)
    )

    # Plot the perfomance coefficients vs advance ratio,
    # for different collective pitches.
    for i, _ in enumerate(theta_75_range):
        axs["eta"].plot(
            J[i], eta[i],
            marker=".", linewidth=0.5, label=f"{theta_75deg_range[i]} (°)"
        )
        axs["thrust"].plot(J[i], C_T[i], marker=".", linewidth=0.5)
        axs["power"].plot(J[i], C_P[i], marker=".", linewidth=0.5)

    axs["eta"].set_xlabel("Advance ratio")
    axs["eta"].set_ylabel("Propulsive efficiency")
    axs["thrust"].set_xlabel("Advance ratio")
    axs["thrust"].set_ylabel("Thrust coefficient")
    axs["power"].set_xlabel("Advance ratio")
    axs["power"].set_ylabel("Power coefficient")

    fig.show()
