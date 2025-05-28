import propu.constant as cst
from propu import jet
from propu.iteralg import IterTable


def main():
    print("Solving exo 4.2: Nozzle change")
    # NOTE: see complement (8)

    ## Problem data and preliminaries

    # Statement data
    M_0 = 0.75
    h = 0  # [m]
    mdot = cst.mconv(165, "lb/s", "kg/s")
    pi_c = 15
    eta_s_c = 0.88
    lhv = cst.mconv(17_800, "BTU/lb", "J/kg")
    T0_3 = cst.mconv(2500, "degR", "K")
    eta_cc = 0.98
    pi_cc = 0.95
    eta_s_t = 0.915
    rr = 0.92
    eta_m = 0.995

    # ISA/SLS assumptions
    sl_state = cst.get_isa(h)
    p_0 = sl_state.p
    T_0 = sl_state.T
    a_0 = sl_state.a
    v_0 = M_0 * a_0
    cp_0 = cst.lerp_cp(T_0, 0)
    g_0 = cp_0 / (cp_0 - cst.R_air)
    T0_0 = T_0 * jet.T_static2total(M_0, g_0)
    p0_0 = p_0 * jet.p_static2total(M_0, g_0)

    # Instantiate the iteration tables
    table_12 = IterTable(("Iter", "cp_12", "T0_2", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_23 = IterTable(("Iter", "cp_3r", "mdot_f", "far"), ("", "(J/(kg*K))", "(kg/s)", ""))
    table_34 = IterTable(("Iter", "cp_34", "T0_4", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_conv = IterTable(
        ("Iter", "cp_5s5t", "Ts_5", "gamma", "NPRC"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    table_div = IterTable(
        ("Iter", "cp_div", "T_5_div", "gamma", "M_5"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    n_iter = 4  # Number of iteration used in iterative procedures

    ## 1. Fuel mass flow rate

    T0_1 = T0_0  # Adiabatic ram compression

    # Find T0_2, cp_12, g_12 iteratively, by using
    # the compressor isentropic efficiency equation.
    g_12 = cst.gamma_air  # initial guess
    cp_12, T0_2, g_12 = jet.compressor_s(T0_1, pi_c, eta_s_c, g_12, table_12)

    # Find mdot_f and far iteratively, by using
    # the combustion equation.
    far = 0.02  # Initial guess
    mdot_f, far = jet.combustion_chamber(T0_2, T0_3, eta_cc, lhv, mdot, far, table_23, n_iter)

    # Mass flow rate at the burner exit [kg/s]
    mdot_b = mdot + mdot_f

    ## 2. Flow condition after turbine nozzle (station 4)

    # Compressor and turbine power via global power balance
    P_c = mdot * cp_12 * (T0_2 - T0_1)
    P_t = P_c / eta_m

    # Find T0_4, cp_34, g_34 iteratively
    cp_34 = cst.lerp_cp(T0_3, far)  # cp_34 guess
    cp_34, T0_4, g_34 = jet.turbine(T0_3, P_t, mdot, far, cp_34, table_34)

    # Find p0_4 from turbine isentropic condition
    p0_3 = p0_0 * rr * pi_c * pi_cc  # chained compressions
    pi_t = (1 - (1 - T0_4 / T0_3) / eta_s_t) ** (-g_34 / (g_34 - 1))
    p0_4 = p0_3 / pi_t

    ## 3. Nozzle comparision

    T0_5 = T0_4  # Adiabatic expansion w/o work
    p0_5 = p0_4  # Isentropic expansion w/o work
    npr = p0_4 / p_0  # Nozzle pressure ratio

    ## 3.1 Converging nozzle

    # NOTE:
    # The NPRC value always lands approximately in the range [1.83 − 1.89].
    # if npr < 1.83:
    #    no need to iterate to find sonic conditions, it will be a waste of time.
    #    The nozzle is adapted.
    # if 1.83 < npr < 1.89:
    #    Need to iterate on sonic condition, to find a more precise value of nprc.
    #    Then see if actually chocked or not.
    # if npr > 1.89:
    #    Nozzle is chocked for sure, but in that case we still need to iterate
    #    to find the sonic conditions.

    # Check nozzle pressure ratio.
    # print(f"{npr=:.4g}")
    # => 5.127
    # => Chocked for sure.

    # NOTE:
    # The computed NPR appears to be well over 1.9
    # => it can safely be assumed that the nozzle is chocked.
    # There was actually no need to get a more precise value of the NPRC through iterations.
    # This NPRC was computed mainly for curiosity, to see its values through iterations.
    # Anyways, the iterations are still necessary for a chocked nozzle,
    # because we have to compute the other sonic flow parameters at the chocked exhaust,
    # like the exhaust speed and pressure.

    # NOTE:
    # This only compute the quantities for a fictitious isentropic expansion
    # to the sonic speed (Ms_5 = 1), so that the nprc can be evaluated.
    # If the nozzle then appears to be chocked, these sonic quantities
    # then corresponds to the actual conditions at the nozzle exhaust.
    convergent = jet.nozzle_sonic(T0_5, p0_5, far, table_conv)
    rho_5_conv = convergent.rho
    v_5_conv = convergent.a
    p_5_conv = convergent.p
    g_5_conv = convergent.g
    # Nozzle throat area, given the mass flow that needs to pass
    A = mdot_b / (rho_5_conv * v_5_conv)

    # Thrust and SFC
    T_conv = (p_5_conv - p_0) * A + mdot_b * v_5_conv - mdot * v_0
    SFC_conv = mdot_f / T_conv

    ## 3.2. Conv-div nozzle

    # NOTE:
    # Conv-div nozzle are always assumed to be adapted.

    # Find the nozzle exhaust conditions iteratively, by using the isentropic relations.
    divergent = jet.nozzle_adapted(T0_5, npr, far, g_5_conv, table_div)
    v_5_div = divergent.v

    # Thrust and SFC
    T_div = mdot_b * v_5_div - mdot * v_0  # p_5 = p_0
    SFC_div = mdot_f / T_div

    ## Print results

    print("1. Iteration table: stations 1 <-> 2")
    table_12.print()
    print("2. Iteration table: stations 2 <-> 3")
    table_23.print()
    print("3. Iteration table: stations 3 <-> 4")
    table_34.print()
    print("4. Iteration table: convergent nozzle")
    table_conv.print()
    print("5. Iteration table: conv-div nozzle")
    table_div.print()
    print("5. Nozzle pressure ratio")
    print(f"   NPR = {npr:.4g} >> NPRC  =>  clearly chocked")
    print("6. Convergent nozzle")
    print(f"   Thrust: {1e-3 * T_conv:.4g} kN")
    print(f"   SFC:    {10 * 3600 * SFC_conv:.4g} kg/(hr*daN)")
    print("7. Conv-div nozzle")
    print(f"   Thrust: {1e-3 * T_div:.4g} kN")
    print(f"   SFC:    {10 * 3600 * SFC_div:.4g} kg/(hr*daN)")
    print("A. For the next exercise")
    print(f"   {p0_4=:.6g} Pa")
    print(f"   {T0_4=:.6g} K")
    print(f"   {far=:.6g}")
    print(f"   {mdot=:.6g} kg/s")
    print(f"   {mdot_f=:.6g} kg/s")
