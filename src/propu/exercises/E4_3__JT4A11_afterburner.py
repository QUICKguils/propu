import propu.constant as cst
import propu.turbine as turb
from propu.iteralg import IterTable


def main():
    print("Solving exo 4.3: Afterburner added on JT4A-11 turbojet")

    # NOTE:
    # This is the continuation of exercise 4_2__Nozzle_change.
    #
    # Basically, start from the turbine outlet conditions (T0_4, p0_4),
    # add afterburner, then compute nozzle as previously.
    #
    # Conditions at 5 are trivial: T0_5 is given, and p0_5 = p0_4 * pi_ab.
    # Apart from that, the afteburner just add another fuel rate than need to be computed
    # iteratively with the combusion equation, just like for the combustion chamber.

    ## Problem data and preliminaries

    # Statement data
    M_0 = 0.75
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
    # Afterburner data
    eta_ab = 0.89
    pi_ab = 0.97
    T0_5 = 1778  # [K]

    # ISA/SLS assumptions
    sl_state = cst.get_isa(0)
    p_0 = sl_state.p
    T_0 = sl_state.T
    a_0 = sl_state.a
    v_0 = M_0 * a_0
    cp_0 = cst.lerp_cp(T_0, 0)
    g_0 = cp_0 / (cp_0 - cst.R_air)
    T0_0 = T_0 * turb.T_static2total(M_0, g_0)
    p0_0 = p_0 * turb.p_static2total(M_0, g_0)

    # From the previous exercise (E4_2__Nozzle_change)
    p0_4 = 519516.1  # Turbine outlet total pressure [Pa]
    T0_4 = 1050.7725  # Turbine outlet total temperature [K]
    far = 0.020329  # fuel-to-air ratio
    mdot = 74.842741  # Inlet mass flow rate [kg/s]
    mdot_f = 1.5214  # Combustion chamber fuel rate [kg/s]
    mdot_b = mdot + mdot_f

    # Instantiate the iteration tables
    table_45 = IterTable(("Iter", "cp_5r", "mdot_f_ab", "far"), ("", "(J/(kg*K))", "(kg/s)", ""))
    table_conv = IterTable(
        ("Iter", "cp_conv", "Ts_6", "gamma", "NPRC"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    table_div = IterTable(
        ("Iter", "cp_div", "T_6_div", "gamma", "M_6"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    n_iter = 4  # Number of iteration used in iterative procedures

    ## Fuel rate in the afterburner

    # Find mdot_f_ab, far_ab, and cp_5r iteratively
    cp_4r = cst.lerp_cp((T0_4 + cst.T_ref) / 2, far)
    far_ab = far  # far guess
    mdot_f_ab = far_ab * mdot - mdot_f  # resulting mdot_f_ab guess
    for iter in range(n_iter):
        cp_5r = cst.lerp_cp((T0_5 + cst.T_ref) / 2, far_ab)
        table_45.add_row(iter, cp_5r, mdot_f_ab, far_ab)  # keep track of iterations
        # Update estimates with combution equation.
        # I could've take time to isolate mdot_f_ab on one side, but this converges anyways.
        mdot_f_ab = (
            (mdot_b + mdot_f_ab) * cp_5r * (T0_5 - cst.T_ref) - mdot_b * cp_4r * (T0_4 - cst.T_ref)
        ) / (eta_ab * lhv)
        far_ab = (mdot_f + mdot_f_ab) / mdot

    ## 3. Nozzle comparision

    # NOTE:
    # This is quite a copy paste of previous exercise.
    # Extended explaination are provided there.

    p0_5 = p0_4 * pi_ab  # Small afterburner pressure loss

    T0_6 = T0_5  # Adiabatic expansion w/o work
    p0_6 = p0_5  # Isentropic expansion w/o work
    npr = p0_5 / p_0  # Nozzle pressure ratio

    ## 3.1 Converging nozzle

    # Check nozzle pressure ratio
    # print(f"{npr=:.4g}")
    # => 4.957
    # => Chocked for sure.

    # Find sonic conditions iteratively
    Ms_6 = 1  # By def. of sonic conditions
    cp_conv = cst.lerp_cp(T0_5, far_ab)  # cp guess
    for iter in range(n_iter):
        g_conv = cp_conv / (cp_conv - cst.R_air)
        nprc = turb.nprc(g_conv)  # Not mandatory. Computed for curiosity
        Ts_6 = T0_6 / turb.T_static2total(M=Ms_6, g=g_conv)
        table_conv.add_row(iter, cp_conv, Ts_6, g_conv, nprc)  # keep track of iterations
        cp_conv = cst.lerp_cp((Ts_6 + T0_6) / 2, far_ab)  # iteration update

    # Exhaust speed and pressure
    a_6 = (g_conv * cst.R_air * Ts_6) ** 0.5
    v_6 = Ms_6 * a_6  # = a_6, as sonic conditions (Ms_6 = 1)
    p_6 = p0_6 / nprc

    # Area of the nozzle, given the mdot_b that needs to pass
    rho_6 = p_6 / (cst.R_air * Ts_6)
    A = mdot_b / (rho_6 * v_6)

    # Thrust and SFC
    T_conv = (p_6 - p_0) * A + (mdot_b + mdot_f_ab) * v_6 - mdot * v_0
    SFC_conv = (mdot_f + mdot_f_ab) / T_conv

    ## 3.2. Conv-div nozzle

    # NOTE:
    # Conv-div nozzle are always assumed to be adapted.
    # So p_6 = p_0, and npr := p0_6/p_0 = p0_6/p_6.

    g_div = g_conv  # gamma guess
    for iter in range(n_iter):
        T_6_div = T0_6 * npr ** ((1 - g_div) / g_div)  # Isentropic relations
        cp_div = cst.lerp_cp((T0_6 + T_6_div) / 2, far_ab)
        M_6 = (((T0_6 / T_6_div) - 1) * 2 / (g_div - 1)) ** 0.5
        table_div.add_row(iter, cp_div, T_6_div, g_div, M_6)  # keep track of iterations
        g_div = cp_div / (cp_div - cst.R_air)  # iteration update

    # Exhaust speed
    a_6 = (g_div * cst.R_air * T_6_div) ** 0.5
    v_6 = M_6 * a_6

    # Thrust and SFC
    T_div = (mdot_b + mdot_f_ab) * v_6 - mdot * v_0  # p_5 = p_0
    SFC_div = (mdot_f + mdot_f_ab) / T_div

    ## Print results

    print("Iteration table: afterburner")
    table_45.print()
    print("Iteration table: convergent nozzle")
    table_conv.print()
    print("Iteration table: conv-div nozzle")
    table_div.print()
    print("Convergent nozzle")
    print(f"  Thrust: {1e-3 * T_conv:.4g} kN")
    print(f"  SFC:    {10 * 3600 * SFC_conv:.4g} kg/(hr*daN)")
    print("Conv-div nozzle")
    print(f"  Thrust: {1e-3 * T_div:.4g} kN")
    print(f"  SFC:    {10 * 3600 * SFC_div:.4g} kg/(hr*daN)")
