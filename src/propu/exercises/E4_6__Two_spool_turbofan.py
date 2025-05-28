import propu.constant as cst
from propu import jet
from propu.iteralg import IterTable


def main():
    print("Solving exo 4.6: Thrust of a two-spool turbofan")
    # NOTE: see schematic on complement (9)

    ## Problem data and preliminaries

    # Statement data
    M_0 = 0.3
    h = cst.mconv(400, "ft", "m")
    mdot_p = 90  # [kg/s]
    bpr = 6
    pi_f = 2.2
    opr = 36
    T0_4 = 1680  # [K]
    pi_cc = 0.95
    eta_cc = 0.99
    lhv = cst.mconv(41.400, "MJ/kg", "J/kg")
    eta_s_f = 0.92
    eta_s_c = 0.9
    eta_s_t = 0.91
    rr = 0.97

    # ISA/SLS assumptions
    upstream_state = cst.get_isa(h)
    p_0 = upstream_state.p
    T_0 = upstream_state.T
    a_0 = upstream_state.a
    v_0 = M_0 * a_0
    cp_0 = cst.lerp_cp(T_0, 0)
    g_0 = cp_0 / (cp_0 - cst.R_air)
    T0_0 = T_0 * jet.T_static2total(M_0, g_0)
    p0_0 = p_0 * jet.p_static2total(M_0, g_0)

    # Immediatly computable
    mdot_s = bpr * mdot_p  # Secondary flow
    pi_c = opr / pi_f  # HPC pressure ratio

    # Instantiate the iteration tables
    table_12 = IterTable(("Iter", "cp_12", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_23 = IterTable(("Iter", "cp_23", "T0_3", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_34 = IterTable(("Iter", "cp_4r", "mdot_f", "far"), ("", "(J/(kg*K))", "(kg/s)", ""))
    table_46 = IterTable(("Iter", "cp_46", "T0_6", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_7 = IterTable(("Iter", "cp_7", "T_7", "gamma", "M_7"), ("", "(J/(kg*K))", "(K)", "", ""))
    table_8 = IterTable(("Iter", "cp_8", "T_8", "gamma", "NPRC"), ("", "(J/(kg*K))", "(K)", "", ""))

    ## 1. Fuel mass flow rate
    #
    # Determine succesively the total temperatures, via the compressor eta_s equation.
    # Then use the combustion equation to find mdot_f.

    T0_1 = T0_0  # Adiabatic ram compression

    # Find cp_12, T0_2, g_12 iteratively, by using the compressor eta-s equation.
    g_12 = cst.gamma_air  # initial guess
    cp_12, T0_2, g_12 = jet.compressor_s(T0_1, pi_f, eta_s_f, g_12, table_12)

    # Find cp_23, T0_3, g_23 iteratively, by using the compressor eta_s equation.
    g_23 = g_12  # initial guess
    cp_23, T0_3, g_23 = jet.compressor_s(T0_2, pi_c, eta_s_c, g_23, table_23)

    # Find mdot_f and far iteratively, by using the combustion equation.
    far = 0.02  # Initial guess
    mdot_f, far = jet.combustion_chamber(T0_3, T0_4, eta_cc, lhv, mdot_p, far, table_34)

    ## 2. Flow conditions at inlet of primary nozzle
    #
    # Determine T0_6 via turbine power balance.
    # Then use the turbine eta_s equation to find p0_6.

    # Power distribution in the turbofan
    P_f = mdot_p * (1 + bpr) * cp_12 * (T0_2 - T0_1)
    P_c = mdot_p * cp_23 * (T0_3 - T0_2)
    P_t = P_f + P_c

    # Find cp_46, T0_6, g_46 iteratively, by using the turbine power equation.
    cp_46 = cst.lerp_cp(T0_4, far)  # initial guess
    cp_46, T0_6, g_46 = jet.turbine(T0_4, P_t, mdot_p, far, cp_46, table_46)

    p0_4 = p0_0 * rr * pi_f * pi_c * pi_cc  # chained compressions
    pi_t = (1 - (1 - T0_6 / T0_4) / eta_s_t) ** (-g_46 / (g_46 - 1))
    p0_6 = p0_4 / pi_t

    ## 3. Nozzle analyses

    ## 3.1 Primary nozzle

    T0_7 = T0_6  # Adiabatic expansion w/o work
    p0_7 = p0_6  # Isentropic expansion w/o work
    npr_p = p0_6 / p_0  # nozzle pressure ratio

    # Check nozzle pressure ratio.
    # print(f"{npr_p=:.4g}")
    # => 1.162
    # => Adapted, for sure.
    # => No need to compute the sonic conditions.

    # Find cp, T, M at primary exhaust iteratively, by using the isentropic relations.
    g_7 = g_46  # initial guess
    M_7 = T_7 = 0  # TBD
    for iter in range(4):
        T_7 = T0_7 * npr_p ** ((1 - g_7) / g_7)  # Isentropic relations
        cp_7 = cst.lerp_cp((T0_7 + T_7) / 2, far)
        M_7 = (((T0_7 / T_7) - 1) * 2 / (g_7 - 1)) ** 0.5
        table_7.add_row(iter, cp_7, T_7, g_7, M_7)  # keep track of iterations
        g_7 = cp_7 / (cp_7 - cst.R_air)  # iteration update

    # Exhaust speed
    a_7 = (g_7 * cst.R_air * T_7) ** 0.5
    v_7 = M_7 * a_7

    ## 3.2 Secondary nozzle

    T0_8 = T0_2  # Adiabatic expansion w/o work
    p0_2 = p0_0 * rr * pi_f  # chained compressions
    p0_8 = p0_2  # Isentropic expansion w/o work
    npr_s = p0_8 / p_0  # nozzle pressure ratio

    # Check nozzle pressure ratio.
    print(f"{npr_s=:.4g}")
    # => 2.272
    # => Chocked, for sure.
    # => Compute the sonic conditions.

    M_8 = 1  # By def. of sonic conditions
    cp_8 = cst.lerp_cp(T0_2, far)  # initial guess
    nprc = T_8 = g_8 = 0  # TBD
    for iter in range(4):
        g_8 = cp_8 / (cp_8 - cst.R_air)
        nprc = jet.get_nprc(g_8)
        T_8 = T0_8 / jet.T_static2total(M=M_8, g=g_8)
        table_8.add_row(iter, cp_8, T_8, g_8, nprc)  # keep track of iterations
        cp_8 = cst.lerp_cp((T_8 + T0_8) / 2, far)  # iteration update

    # Exhaust speed and pressure
    a_8 = (g_8 * cst.R_air * T_8) ** 0.5
    v_8 = M_8 * a_8  # = a_8, as sonic conditions (Ms_5 = 1)
    p_8 = p0_8 / nprc

    # Area of the nozzle, given the mdot_b that needs to pass
    rho_8 = p_8 / (cst.R_air * T_8)
    A_8 = mdot_p * bpr / (rho_8 * v_8)

    ## 3.3 Thrust and SFC

    T = (mdot_p + mdot_f) * v_7 + mdot_p * bpr * v_8 - mdot_p * (1 + bpr) * v_0 + (p_8 - p_0) * A_8
    SFC = mdot_f / T

    ## Print results

    print("Iteration table: fan")
    table_12.print()
    print("Iteration table: high-pressure compressor")
    table_23.print()
    print("Iteration table: combustion chamber")
    table_34.print()
    print("Iteration table: turbine")
    table_46.print()
    print("Iteration table: primary exhaust")
    table_7.print()
    print("Iteration table: secondary exhaust")
    table_8.print()
    print(f"Thrust: {1e-3 * T:.4g} kN")
    print(f"SFC:    {10 * 3600 * SFC:.4g} kg/(hr*daN)")
