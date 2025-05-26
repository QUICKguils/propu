import math

import propu.constant as cst
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.3: Turbine polytropic efficiency")
    # NOTE: see complement (5)

    ## Problem data and preliminaries

    # Statement data
    Pm = 30 * cst.uconv("MW", "W")
    T0_out = cst.mconv(1473, "degR", "K")
    mdot_a = 164 * cst.uconv("lb/s", "kg/s")
    far = 0.015
    eta_s = 0.93

    # Instantiate an iteration table
    variables = ("Iter", "pi_t", "eta_p", "gamma")
    units = ("", "", "(%)", "")
    table = IterTable(variables, units)

    # Exhaust mass flow rate [kg/s]
    mdot_b = mdot_a * (1 + far)

    ## Resolution

    # Initial guess for g, cp and T0_in
    g = cst.gamma_air
    cp = cst.R_air * g / (g - 1)
    T0_in = T0_out + Pm / (mdot_b * cp)

    for iter in range(5):
        # Compute the resulting quantities,
        # that are function of cp and/or T0_in
        g = cp / (cp - cst.R_air)
        pi_t = (1 - (1 - T0_out / T0_in) / eta_s) ** (-g / (g - 1))
        eta_p = g / (g - 1) * math.log(T0_in / T0_out) / math.log(pi_t)

        # Add the computed data for printing
        table.add_row(iter, pi_t, 100 * eta_p, g)

        # Iteration update
        cp = cst.lerp_cp((T0_in + T0_out) / 2, far)
        T0_in = T0_out + Pm / (mdot_b * cp)

    table.print()
