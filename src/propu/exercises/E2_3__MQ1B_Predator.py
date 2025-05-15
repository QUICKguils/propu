import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize

from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.mplrc import REPORT_TW


def main():
    print("Solving exo 2.3: MQ-1B Predator")

    # 0. Statement data

    Pm = 115 * uconv("hp", "W")  # Mechanical power (Pm = Ps here)
    rho = get_state(0).rho
    D = 6 * uconv("ft", "m")

    # 1. Prepeller efficiency vs. airspeed
    # NOTE: See complement (2)

    Ap_single = np.pi * D**2 / 4

    def f(eta, v_inf, Ap):
        return eta**3 + (eta - 1) * v_inf**3 * 2 * rho * Ap / Pm

    guess = 0.5
    v_range = np.linspace(2, 100, 50)
    eta_range = np.array([optimize.root(f, guess, args=(v, Ap_single)).x for v in v_range])

    # 2. Static thrust vs. propeller size
    # NOTE: in the statement, he means "ft" and not "in"

    D_range = np.linspace(0, 12 * uconv("ft", "m"), 50)
    # D_range = np.linspace(0, 4, 50)  # As in the exo solution
    Ap_range = np.pi * D_range**2 / 4
    Ts_range = (2 * rho * Ap_range * Pm**2) ** (1 / 3)  # Just Eq.(2.36) p.13-2

    # 3. Thrust vs. airspeed
    # NOTE: in the statement, he means "ft" and not "in"

    D_selected = np.array([3, 6, 9, 12]) * uconv("ft", "m")
    Ap_selected = np.pi * D_selected**2 / 4
    thrust_selected = np.zeros((Ap_selected.size, v_range.size))
    for i, Ap in enumerate(Ap_selected):
        eta_selected = np.ravel([optimize.root(f, guess, args=(v, Ap)).x for v in v_range])
        thrust_selected[i, :] = eta_selected * Pm / v_range

    # A. Plots

    fig, axs = plt.subplots(1, 3, figsize=(2 * REPORT_TW, 0.6* REPORT_TW))
    fig.suptitle(f"MQ-1B Predator @ {Pm / 1e3:.2f} kW")

    axs[0].plot(v_range, eta_range)
    axs[0].set_xlabel(r"$v_{inf}$/(m/s)")
    axs[0].set_ylabel(r"$\eta_p$")

    axs[1].plot(D_range, Ts_range)
    axs[1].set_xlabel(r"$D$/m")
    axs[1].set_ylabel(r"$T_s$/N")

    for i, _ in enumerate(Ap_selected):
        axs[2].plot(v_range, thrust_selected[i, :])
    axs[2].set_xlabel(r"$v_{inf}$/(m/s)")
    axs[2].set_ylabel(r"$T$/N")

    fig.show()
