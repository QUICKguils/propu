import math

import propu.constant as c
from propu.iteralg import IterTable


def main():
    print("Solving exo 3.3: Turbine polytropic efficiency")
    # NOTE: see complement (5)

    # Statement data
    Pm = 30 * c.uconv("MW", "W")
    T0_out = c.mconv(1473, "degR", "K")
    mdot_a = 164 * c.uconv("lb/s", "kg/s")
    far = 0.015
    eta_s = 0.93

    # Instantiate a table printer
    variables = ("Iter", "pi_t", "eta_p", "gamma")
    units = ("", "", "(%)", "")
    table = IterTable(variables, units)

    # Mass flow of fuel + air [kg/s]
    mdot_b = mdot_a * (1 + far)

    # Initial guess for g, cp and T0_in
    g = c.gamma_air
    cp = c.R_air * g / (g - 1)
    T0_in = T0_out + Pm / (mdot_b * cp)
    n_iter = 5
    for iter in range(1, n_iter + 1):
        # Compute the resulting quantities,
        # that are function of cp and/or T0_in
        g = cp / (cp - c.R_air)
        pi_t = (1 - (1 - T0_out / T0_in) / eta_s) ** (-g / (g - 1))
        eta_p = g / (g - 1) * math.log(T0_in / T0_out) / math.log(pi_t)

        # Add the computed data for printing
        table.add_row(iter, pi_t, 100 * eta_p, g)

        # Iteration update
        cp = c.lerp_cp((T0_in + T0_out) / 2, far)
        T0_in = T0_out + Pm / (mdot_b * cp)

    table.print()
