import propu.constant as cst
from propu import jet
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.4: Power losses")
    # NOTE: see complement (6)

    ## Problem data and preliminaries

    # Statement data
    pi_c = 35
    mdot_b = cst.mconv(351, "lb/s", "kg/s")
    T0_4 = cst.mconv(725, "degF", "K")
    P_s = cst.mconv(51.2, "MW", "W")
    HR = cst.mconv(8210, "BTU/kWh", "")
    lhv = cst.mconv(47_110, "kJ/kg", "J/kg")
    eta_s_c = 0.88
    eta_cc = 1

    # Mass flow rates [kg/s]
    mdot_f = HR * P_s / (eta_cc * lhv)
    mdot_a = mdot_b - mdot_f
    far = mdot_f / mdot_a

    # Instantiate the iteration tables
    table_12 = IterTable(("Iter", "cp_c", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_13 = IterTable(("Iter", "cp_13", "T0_3"), ("", "(J/(kg*K))", "(K)"))

    ## 1. Stations 1 -> 2 : find T0_2 iteratively

    # Find cp_12, T0_2, g_12 iteratively, by using the compressor eta-s equation.
    g_12 = cst.gamma_air  # initial guess
    cp_12, T0_2, g_12 = jet.compressor_s(cst.T_ref, pi_c, eta_s_c, g_12, table_12)

    # 2. Stations 1 -> 3 : find T0_3 iteratively

    g_13 = cst.gamma_air  # gamma guess
    cp_13 = cst.R_air * g_13 / (g_13 - 1)  # cp guess
    for iter in range(5):
        T0_3 = cst.T_ref + (eta_cc * lhv * mdot_f + mdot_a * cp_12 * (T0_2 - cst.T_ref)) / (
            mdot_b * cp_13
        )
        table_13.add_row(iter, cp_13, T0_3)  # keep track of iterations
        cp_13 = cst.lerp_cp((T0_3 + cst.T_ref) / 2, far)  # iteration update

    ## 3. Powers from global energy conservation

    P_c = mdot_a * cp_12 * (T0_2 - cst.T_ref)
    cp_t = cst.lerp_cp((T0_3 + T0_4) / 2, far)
    P_t = mdot_b * cp_t * (T0_3 - T0_4)
    P_fa = P_t - P_c - P_s

    ## Print results

    print("1. Iteration table for T0_2")
    table_12.print()
    print("2. Iteration table for T0_3")
    table_13.print()
    print(
        "3. Gobal power results",
        f"  Shaft:      P_s  = {1e-6 * P_s:.4g} MW (given)",
        f"  Compressor: P_c  = {1e-6 * P_c:.4g} MW",
        f"  Turbine:    P_t  = {1e-6 * P_t:.4g} MW",
        f"  Loss:       P_fa = {1e-6 * P_fa:.4g} MW",
        sep="\n",
    )
