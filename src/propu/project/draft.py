import numpy as np

from propu.project import statement as stm
from propu.project import part1
from propu import bemt


def plot_pitches() -> None:
    import matplotlib.pyplot as plt
    from propu.mplrc import REPORT_TW
    fig, axs = plt.subplots(2, 2, figsize=(REPORT_TW, REPORT_TW))
    oper = part1.get_operating_conditions()
    aoa_opt = find_optimal_aoa()
    aoa_opt_deg = np.rad2deg(aoa_opt)
    print(f"{aoa_opt_deg=:.2f}")

    def plot(prop, ax):
        ax.plot(
            prop.geometry.stations/prop.geometry.span,
            np.rad2deg(prop.geometry.pitches),
            label="Blade pitch"
        )
        ax.plot(
            sol.r_dist/prop.geometry.span,
            np.rad2deg(np.pi/2 + aoa_opt + sol.beta_dist),
            label="Optimal pitch"
        )
        ax.set_title(f"{prop.pretty_name}")
        ax.set_xlabel(r"$r/R$")
        ax.set_ylabel(r"$\bar{\chi}$")

    for i, prop_keyword in enumerate(stm.APC_dict):
        prop = stm.create_apc_propeller(prop_keyword)
        sol = bemt.bem(prop, oper)
        plot(prop, axs.flat[i])

    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='outside upper center', ncols=2)
    fig.show()


def find_optimal_aoa() -> float:
    arg_aoa = int(np.mean(np.argmax(stm.airfoil.cl / stm.airfoil.cd, axis=0)))
    # arg_aoa = np.argmin(stm.airfoil.cd, axis=0)[0]
    return stm.airfoil.aoa[arg_aoa]


if __name__ == '__main__':
    plot_pitches()
