import propu.constant as cst
from propu import jet
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.5: Jet fuel starter")
    # NOTE: see complement (7)

    ## Problem data and preliminaries

    # Statement data
    P_s = cst.mconv(172, "kW", "W")
    sfc = cst.mconv(878, "g/kWh", "kg/J")
    mdot_a = 1.75  # [kg/s]
    T0_3 = cst.mconv(983, "degC", "K")
    eta_s_c = 0.85  # Aww badd :(
    pi_c = 2.86
    lhv = cst.mconv(42_979, "kJ/kg", "J/kg")

    # ISA/SLS assumptions
    T0_1 = cst.T_ref

    # Mass flow rates
    mdot_f = sfc * P_s
    mdot_b = mdot_a + mdot_f
    far = mdot_f / mdot_a

    # Instantiate the iteration tables
    table_12 = IterTable(("Iter", "cp_12", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_34 = IterTable(("Iter", "cp_34", "T0_4"), ("", "(J/(kg*K))", "(K)"))

    ## 1. Stations 1 -> 2 : find T0_2 iteratively
    #
    # A bit overkill to iterate here,
    # as the temperature gap is expected to be not that high,
    # for a compressor with quite low pressure ratio.

    # Find cp_12, T0_2, g_12 iteratively, by using the compressor eta-s equation.
    g_12 = cst.gamma_air  # initial guess
    cp_12, T0_2, g_12 = jet.compressor_s(T0_1, pi_c, eta_s_c, g_12, table_12)

    ## 2. Between the two turbines : find T0_4 iteratively

    # Turbine power from global power balances
    P_c = mdot_a * cp_12 * (T0_2 - T0_1)
    P_t = P_c  # no P_s, as not already crossed the second turbine

    # Find turbine outlet conditions iteratively, by using the turbine power equation.
    _ = jet.turbine(T0_3, P_t, mdot_a, far, cp_12, table_34)  # just side effect of completing the table

    ## 3. Efficiency of the combustion chamber

    # Just apply definition of eta_cc, from combustion equation
    cp_3r = cst.lerp_cp((T0_3 + T0_1) / 2, far)
    cp_2r = cst.lerp_cp((T0_2 + T0_1) / 2, 0)
    eta_cc = (mdot_b * cp_3r * (T0_3 - T0_1) - mdot_a * cp_2r * (T0_2 - T0_1)) / (mdot_f * lhv)

    ## Print results

    print("1. Iteration table for T0_2")
    table_12.print()
    print("2. Iteration table for T0_4")
    table_34.print()
    print(f"3. Combustion chamber efficiency: {100 * eta_cc:.4g} %")
