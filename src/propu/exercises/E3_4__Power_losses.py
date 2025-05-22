import propu.constant as c
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.4: Power losses")
    # NOTE: see complement (6)

    # Statement data
    pi_c = 35
    mdot_b = c.mconv(351, "lb/s", "kg/s")
    T0_4 = c.mconv(725, "degF", "K")
    P_s = c.mconv(51.2, "MW", "W")
    HR = c.mconv(8210, "BTU/kWh", "")
    lhv = c.mconv(47_110, "kJ/kg", "J/kg")
    eta_s_c = 0.88
    eta_cc = 1

    # Mass flow rates [kg/s]
    mdot_f = HR * P_s / (eta_cc * lhv)
    mdot_a = mdot_b - mdot_f
    far = mdot_f / mdot_a

    # 1. Stations 1 -> 2 : find T0_2 iteratively

    # Instantiate a table printer
    table_12 = IterTable(("Iter", "cp_c", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))

    g_12 = c.gamma_air
    n_iter = 5
    for iter in range(1, n_iter + 1):
        T0_2 = c.T_ref * (1 + 1 / eta_s_c * (pi_c ** ((g_12 - 1) / g_12) - 1))
        cp_12 = c.lerp_cp((T0_2 + c.T_ref) / 2, 0)

        # Add the computed data for printing
        table_12.add_row(iter, cp_12, T0_2, g_12)

        # Iteration update
        g_12 = cp_12 / (cp_12 - c.R_air)

    print("Iteration table for T0_2")
    table_12.print()

    # 2. Stations 1 -> 3 : find T0_3 iteratively

    # Instantiate a table printer
    table_13 = IterTable(("Iter", "cp_13", "T0_3"), ("", "(J/(kg*K))", "(K)"))

    g_13 = c.gamma_air
    cp_13 = c.R_air * g_13 / (g_13 - 1)
    n_iter = 5
    for iter in range(1, n_iter + 1):
        T0_3 = c.T_ref + (eta_cc * lhv * mdot_f + mdot_a * cp_12 * (T0_2 - c.T_ref)) / (
            mdot_b * cp_13
        )

        # Add the computed data for printing
        table_13.add_row(iter, cp_13, T0_3)

        # Iteration update
        cp_13 = c.lerp_cp((T0_3 + c.T_ref) / 2, far)

    print("Iteration table for T0_3")
    table_13.print()

    # 3. Powers from global energy conservation
    P_c = mdot_a * cp_12 * (T0_2 - c.T_ref)
    cp_t = c.lerp_cp((T0_3 + T0_4) / 2, far)
    P_t = mdot_b * cp_t * (T0_3 - T0_4)
    P_fa = P_t - P_c - P_s

    print(
        "Gobal power results:",
        f"  Shaft:      P_s  = {1e-6 * P_s:.4g} MW (given)",
        f"  Compressor: P_c  = {1e-6 * P_c:.4g} MW",
        f"  Turbine:    P_t  = {1e-6 * P_t:.4g} MW",
        f"  Loss:       P_fa = {1e-6 * P_fa:.4g} MW",
        sep="\n",
    )
