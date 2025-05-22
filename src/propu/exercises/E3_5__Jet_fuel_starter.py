import propu.constant as c
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.5: Jet fuel starter")

    ## Problem data and preliminaries {{{1

    # Statement data
    P_s = c.mconv(172, "kW", "W")
    sfc = c.mconv(878, "g/kWh", "kg/J")
    mdot_a = 1.75  # [kg/s]
    T0_3 = c.mconv(983, "degC", "K")
    eta_s_c = 0.85  # Aww badd :(
    pi_c = 2.86
    lhv = c.mconv(42_979, "kJ/kg", "J/kg")

    # ISA/SLS assumptions
    T0_1 = c.T_ref
    p0_1 = c.p_ref

    # Mass flow rates
    mdot_f = sfc * P_s
    mdot_b = mdot_a + mdot_f
    far = mdot_f / mdot_a

    # Instantiate the iteration tables
    table_12 = IterTable(("Iter", "cp_12", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_34 = IterTable(("Iter", "cp_34", "T0_4"), ("", "(J/(kg*K))", "(K)"))

    ## 1. Stations 1 -> 2 : find T0_2 iteratively {{{1
    #
    # A bit overkill to iterate here,
    # as the temperature gap is expected to be not that high,
    # for a compressor with quite low pressure ratio.

    g_12 = c.gamma_air  # guesstimate for gamma
    for iter in range(3):
        T0_2 = c.T_ref * (1 + 1 / eta_s_c * (pi_c ** ((g_12 - 1) / g_12) - 1))
        cp_12 = c.lerp_cp((T0_2 + c.T_ref) / 2, 0)
        table_12.add_row(iter, cp_12, T0_2, g_12)  # keep track of interations
        g_12 = cp_12 / (cp_12 - c.R_air)  # iteration update

    ## 2. Find T0_4 between the two turbines {{{1

    # Turbine power from global power balances
    P_c = mdot_a * cp_12 * (T0_2 - T0_1)
    P_t = P_c  # no P_s, as not already crossed the second turbine

    cp_34 = c.R_air * c.gamma_air / (c.gamma_air - 1)  # guesstimates for cp
    for iter in range(3):
        T0_4 = T0_3 - P_t / (mdot_b * cp_34)
        table_34.add_row(iter, cp_34, T0_4)  # keep track of interations
        cp_34 = c.lerp_cp((T0_3 + T0_4) / 2, far)  # iteration update

    ## 3. Efficiency of the combustion chamber {{{1

    # Just apply definition of eta_cc, from combustion equation
    cp_3r = c.lerp_cp((T0_3 + T0_1) / 2, far)
    cp_2r = c.lerp_cp((T0_2 + T0_1) / 2, 0)
    eta_cc = (mdot_b * cp_3r * (T0_3 - T0_1) - mdot_a * cp_2r * (T0_2 - T0_1)) / (mdot_f * lhv)

    ## A. Print results {{{1

    print("1. Iteration table for T0_2")
    table_12.print()
    print("2. Iteration table for T0_4")
    table_34.print()
    print(f"3. Combustion chamber efficiency: {100 * eta_cc:.4g} %")
