import numpy as np

import propu.constant as cst
from propu.iteralg import IterTable


def main():
    print("Solving exo 4.1: JT4A-11 turbojet")

    ## Problem data and preliminaries

    # Statement data
    altitude = cst.mconv(25_000, "ft", "m")
    lhv = cst.mconv(42_979, "kJ/kg", "J/kg")
    M_list = np.arange(0.0, 0.9, 0.1)

    # Estimated from the graph (was a pain to read ngl)
    T_list = (
        np.array([7000, 6800, 6600, 6500, 6500, 6600, 6700, 6800, 7000, 7500, 7900])
        * cst.uconv("lb", "kg")
        * cst.ge
    )
    mdot_list = (
        np.array([105, 106, 107, 112, 115, 120, 125, 135, 147, 156, 170])
        * cst.uconv("lb/s", "kg/s")
    )
    TSFC_list = (
        np.array([0.73, 0.76, 0.8, 0.82, 0.86, 0.87, 0.91, 0.95, 0.97, 1.01, 1.05])
        * cst.uconv("lb/(hr*lb)", "kg/(s*kg)")
        / cst.ge
    )

    # ISA/SLS assumptions
    a_0 = cst.get_isa(altitude).a  # Upstream airspeed [m/s]

    # Instantiate an iteration table
    var_names = ("M_0", "v_j", "eta_h", "eta_p", "eta_g", "TSFC")
    unit_names = ("", "(m/s)", "(%)", "(%)", "(%)", "(kg/(s*N))")
    table = IterTable(var_names, unit_names)

    ## Resolution

    for M, T, mdot, TSFC in zip(M_list, T_list, mdot_list, TSFC_list):
        v_0 = M * a_0
        mdot_f = TSFC * T

        # Just wrangle the performance parameters formulae
        # See (4.47) to (4.51), p.38-1
        v_j = T / mdot + v_0  # Pressure neglected, as isentropic expansion is assumed
        eta_t = mdot * (v_j**2 - v_0**2) / (2 * mdot_f * lhv)
        eta_p = 2 * v_0 / (v_0 + v_j)
        eta_g = eta_t * eta_p
        table.add_row(M, v_j, 100 * eta_t, 100 * eta_p, 100 * eta_g, TSFC)

    table.print()
