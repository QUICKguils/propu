"""Answer the first question of the exam."""

import propu.constant as cst
from propu import jet
from propu.iteralg import IterTable


def main():
    """Execute the first question of the exam."""

    ## Problem data and preliminaries

    # Statement data
    mdot_p = 77  # Air mass flow rate [kg/s]
    mdot_f_ab = 3.2  # Afterburner fuel mass flow rate [kg/s]
    opr = 13.5  # Overall pressure ratio
    T0_4 = 1210  # Turbine inlet temperature [K]
    lhv = cst.mconv(42.8, "MJ/kg", "J/kg")  # Fuel lower heating value
    eta_m = 0.99  # Shaft efficiency
    rr = 0.99  # Ram recovery  factor
    eta_p_c = 0.98  # Compressor polytropic efficiency
    eta_s_t = 0.98  # Turbine isentropic efficiency
    eta_cc = 0.97  # Combustion chamber efficiency
    pi_cc = 0.97  # Combustion chamber pressure ratio
    eta_ab = 0.97  # Afterburner combustion efficiency
    pi_ab = 0.95  # Afterburner pressure ratio

    # ISA/SLS/cruise assumptions
    v_0 = cst.mconv(590.3, "mph", "m/s")
    h = cst.mconv(35_000, "ft", "m")
    # Resulting upstream flow conditions
    upstream_state = cst.get_isa(h)
    p_0 = upstream_state.p
    T_0 = upstream_state.T
    a_0 = upstream_state.a
    M_0 = v_0 / a_0
    cp_0 = cst.lerp_cp(T_0, 0)
    g_0 = cp_0 / (cp_0 - cst.R_air)
    T0_0 = T_0 * jet.T_static2total(M_0, g_0)
    p0_0 = p_0 * jet.p_static2total(M_0, g_0)

    # Instantiate the iteration tables
    table_23 = IterTable(("Iter", "cp_23", "T0_3", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_34 = IterTable(("Iter", "cp_4r", "mdot_f", "far"), ("", "(J/(kg*K))", "(kg/s)", ""))
    table_45 = IterTable(("Iter", "cp_45", "T0_5", "gamma"), ("", "(J/(kg*K))", "(K)", ""))
    table_56 = IterTable(("Iter", "cpr_6", "T0_6", "g_6r"), ("", "(J/(kg*K))", "(K)", ""))
    table_7_wet = IterTable(
        ("Iter", "cp_7", "T_7", "gamma", "NPRC"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    table_8_wet = IterTable(
        ("Iter", "cp_8", "T_8", "gamma", "M_8"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    table_7_dry = IterTable(
        ("Iter", "cp_7", "T_7", "gamma", "NPRC"), ("", "(J/(kg*K))", "(K)", "", "")
    )
    table_8_dry = IterTable(
        ("Iter", "cp_8", "T_8", "gamma", "M_8"), ("", "(J/(kg*K))", "(K)", "", "")
    )

    ## 1. Establish the total conditions

    # 1.1 Inlet
    # Assumed to be an isentropic compression of the air.
    # So the total temperature and pressure are simply conserved.
    T0_1 = T0_0
    p0_1 = p0_0

    # 1.2 Diffuser
    # Assumed to be an adiabatic compression, so that T0 is conserved.
    # But small pressure loss, so that p0 is not totally conserved.
    T0_2 = T0_1
    p0_2 = p0_1 * rr

    # 1.3 Compressor
    p0_3 = p0_2 * opr
    # cp23, T0_3 and g_23 are found iteratively, by using the compressor eta-p equation.
    g_23 = g_0  # initial guess
    cp_23, T0_3, g_23 = jet.compressor_p(T0_2, opr, eta_p_c, g_23, table_23)

    # 1.4 Combustion chamber
    p0_4 = p0_3 * pi_cc
    # T0_4 is known (TIT is a design parameter that is often given).
    # But we still need to find mdot_f and far iteratively, by using the combustion equation.
    far = 0.02  # Initial guess
    mdot_f, far = jet.combustion_chamber(T0_3, T0_4, eta_cc, lhv, mdot_p, far, table_34)
    # The fuel-to-air ratio in the afterburner can now be computed
    far_ab = (mdot_f + mdot_f_ab) / mdot_p

    # 1.5 Turbine
    # First determine the power distribution in the jet
    P_c = mdot_p * cp_23 * (T0_3 - T0_2)
    P_t = P_c / eta_m
    # Then cp45, T0_5 and g_45 are found iteratively, by using the turbine power equation.
    cp_45 = cst.lerp_cp(T0_4, far)  # initial guess
    cp_45, T0_5, g_45 = jet.turbine(T0_4, P_t, mdot_p, far, cp_45, table_45)
    # The pressure ratio of the turbine can then be found, and hence the outlet total pressure.
    pi_t = (1 - (1 - T0_5 / T0_4) / eta_s_t) ** (-g_45 / (g_45 - 1))
    p0_5 = p0_4 / pi_t

    # 1.6 Afterburner
    p0_6 = p0_5 * pi_ab
    # T0_6 is found iteratively, by using the combustion equation.
    cp_6r, T0_6, g_6r = jet.afterburner(T0_5, eta_ab, lhv, mdot_p, mdot_f, mdot_f_ab, table_56)

    ## 2. Nozzle analysis

    # NOTE:
    # In this course, nozzle are always assumed to isentropically expand the gas.
    # Thus, the total temperature and total pressure are conserved throughout the nozzle.
    # Moreover, the conv-div nozzle are always considered to be adapted.
    # This means that the exhaust pressure p_8 is equal to the back-pressure,
    # or atmospheric pressure p_0.
    # Note also that the sonic conditions are reached in the throat (otherwise, the flow
    # has to remain subsonic everywhere in the nozzle, and it would not make sense to use
    # a conv-div nozzle that will decelerate this subsonic flow in the divergent section).

    # 2.1 Wet conditions

    T0_7_wet = T0_8_wet = T0_6  # Adiabatic expansion w/o work
    p0_7_wet = p0_8_wet = p0_6  # Isentropic expansion w/o work
    npr_wet = p0_6 / p_0  # nozzle pressure ratio

    # Find the nozzle throat conditions iteratively, by using the isentropic relations.
    sonic_wet = jet.nozzle_sonic(T0_7_wet, p0_7_wet, far_ab, table_7_wet)
    rho_7_wet = sonic_wet.rho
    v_7_wet = sonic_wet.a
    g_7_wet = sonic_wet.g
    # Nozzle throat area, given the mass flow that needs to pass
    A_7_wet = (mdot_p + mdot_f + mdot_f_ab) / (rho_7_wet * v_7_wet)

    # Find the nozzle exhaust conditions iteratively, by using the isentropic relations.
    exhaust_wet = jet.nozzle_adapted(T0_8_wet, npr_wet, far_ab, g_7_wet, table_8_wet)
    v_8_wet = exhaust_wet.v
    # Thrust (no contribution from pressure difference, as the nozzle is adaped)
    thrust_wet = (mdot_p + mdot_f + mdot_f_ab) * v_8_wet - mdot_p * v_0

    # 2.2 Dry conditions

    # NOTE:
    # We can simply ignore the afterburner, and use the outlet turbine conditions
    # (station 5) as the inlet of the nozzle.
    # One should however be careful to not take into account mdot_f_ab anymore.

    T0_7_dry = T0_8_dry = T0_5  # Adiabatic expansion w/o work
    p0_7_dry = p0_8_dry = p0_5  # Isentropic expansion w/o work
    npr_dry = p0_5 / p_0  # nozzle pressure ratio

    # Find the nozzle throat conditions iteratively, by using the isentropic relations.
    sonic_dry = jet.nozzle_sonic(T0_7_dry, p0_7_dry, far_ab, table_7_dry)
    # Extract quantities useful for the following
    rho_7_dry = sonic_dry.rho
    v_7_dry = sonic_dry.a
    g_7_dry = sonic_dry.g
    # Nozzle throat area, given the mass flow that needs to pass
    A_7_dry = (mdot_p + mdot_f) / (rho_7_dry * v_7_dry)

    # Find the nozzle exhaust conditions iteratively, by using the isentropic relations.
    exhaust_dry = jet.nozzle_adapted(T0_8_dry, npr_dry, far, g_7_dry, table_8_dry)
    v_8_dry = exhaust_dry.v
    # Thrust (no contribution from pressure difference, as the nozzle is adaped)
    thrust_dry = (mdot_p + mdot_f) * v_8_dry - mdot_p * v_0

    ## 3. Performance analysis

    ## 3.1 Wet conditions

    # Powers
    P_t_wet = (mdot_f + mdot_f_ab) * lhv
    P_p_wet = thrust_wet * v_0
    P_m_wet = 0.5 * ((mdot_p + mdot_f + mdot_f_ab) * v_8_wet**2 - mdot_p * v_0**2)
    # Losses
    thermal_loss_wet = P_t_wet - P_m_wet
    kinetic_loss_wet = P_m_wet - P_p_wet
    # Coefficients
    eta_prop_wet = P_p_wet / P_m_wet
    eta_thermal_wet = P_m_wet / P_t_wet
    eta_overall_wet = eta_prop_wet * eta_thermal_wet

    ## 3.2 Dry conditions

    # Powers
    P_t_dry = mdot_f * lhv
    P_p_dry = thrust_dry * v_0
    P_m_dry = 0.5 * ((mdot_p + mdot_f) * v_8_dry**2 - mdot_p * v_0**2)
    # Losses
    thermal_loss_dry = P_t_dry - P_m_dry
    kinetic_loss_dry = P_m_dry - P_p_dry
    # Coefficients
    eta_prop_dry = P_p_dry / P_m_dry
    eta_thermal_dry = P_m_dry / P_t_dry
    eta_overall_dry = eta_prop_dry * eta_thermal_dry

    ## Print results

    print("\n1. Total quantities throughout the jet")

    print("Station 0: Upstream conditions")
    print(f"p0_0 = {1e-3 * p0_0:.4g} kPa")
    print(f"T0_0 = {T0_0:.4g} K")
    print("Station 1: intake")
    print(f"p0_1 = {1e-3 * p0_1:.4g} kPa")
    print(f"T0_1 = {T0_1:.4g} K")
    print("Station 2: compressor inlet")
    print(f"p0_2 = {1e-3 * p0_2:.4g} kPa")
    print(f"T0_2 = {T0_2:.4g} K")
    print("Station 3: compressor outlet")
    print(f"p0_3 = {1e-3 * p0_3:.4g} kPa")
    print(f"T0_3 = {T0_3:.4g} K")
    table_23.print()
    print("Station 4: turbine inlet")
    print(f"p0_4 = {1e-3 * p0_4:.4g} kPa")
    print(f"T0_4 = {T0_4:.4g} K")
    table_34.print()
    print("Station 5: turbine outlet")
    print(f"p0_5 = {1e-3 * p0_5:.4g} kPa")
    print(f"T0_5 = {T0_5:.4g} K")
    table_45.print()
    print("Station 6: Afterburner outlet")
    print(f"p0_6 = {1e-3 * p0_6:.4g} kPa")
    print(f"T0_6 = {T0_6:.4g} K")
    table_56.print()

    print("\n2. Nozzle analysis")

    print("\n2.1 Wet conditions")
    print("2.1.1 Conditions at nozzle throat")
    print(f"Throat area: {A_7_wet:.4g} m**2")
    table_7_wet.print()
    print("2.1.2 Conditions at nozzle exhaust")
    print(f"Exhaust speed: {v_8_wet:.4g} m/s")
    print(f"Thrust: {1e-3 * thrust_wet:.4g} kN")
    table_8_wet.print()

    print("\n2.2 Dry conditions")
    print("2.1.1 Conditions at nozzle throat")
    print(f"Throat area: {A_7_dry:.4g} m**2")
    table_7_dry.print()
    print("2.1.2 Conditions at nozzle exhaust")
    print(f"Exhaust speed: {v_8_dry:.4g} m/s")
    print(f"Thrust: {1e-3 * thrust_dry:.4g} kN")
    table_8_dry.print()

    print("\n3. Performance analysis")

    print("\n3.1 Wet conditions")
    print(f"Thermal power: {1e-6 * P_t_wet:.4g} MW")
    print(f"Propulsive power: {1e-6 * P_p_wet:.4g} MW")
    print(f"Air kinetic power: {1e-6 * P_m_wet:.4g} MW")
    print(f"Thermal loss: {1e-6 * thermal_loss_wet:.4g} MW")
    print(f"Air kinetic loss: {1e-6 * kinetic_loss_wet:.4g} MW")
    print(f"Propulsive efficiency: {eta_prop_wet:.4g}")
    print(f"Thermal efficiency: {eta_thermal_wet:.4g}")
    print(f"Overall efficiency: {eta_overall_wet:.4g}")

    print("\n3.2 Dry conditions")
    print(f"Thermal power: {1e-6 * P_t_dry:.4g} MW")
    print(f"Propulsive power: {1e-6 * P_p_dry:.4g} MW")
    print(f"Air kinetic power: {1e-6 * P_m_dry:.4g} MW")
    print(f"Thermal loss: {1e-6 * thermal_loss_dry:.4g} MW")
    print(f"Air kinetic loss: {1e-6 * kinetic_loss_dry:.4g} MW")
    print(f"Propulsive efficiency: {eta_prop_dry:.4g}")
    print(f"Thermal efficiency: {eta_thermal_dry:.4g}")
    print(f"Overall efficiency: {eta_overall_dry:.4g}")
