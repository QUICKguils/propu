"""Jet engine cycle analyses."""

import propu.constant as cst
from propu.iteralg import IterTable


def _f(M, g):
    """Function f(M), from appendix C."""
    return 1 + (g - 1) / 2 * M**2


def _F(M, g):
    """Function F(M), from appendix C."""
    return g**0.5 * M * _f(M, g) ** (-(g + 1) / (2 * g - 2))


def _G(M, g):
    """Function G(M), from appendix C."""
    return (1 + g * M**2) * _f(M, g) ** (-g / (g - 1))


def T_static2total(M, g):
    """Static to total temperature conversion factor."""
    return _f(M, g)


def p_static2total(M, g):
    """Static to total pressure conversion factor."""
    return _f(M, g) ** (g / (g - 1))


def rho_static2total(M, g):
    """Static to total density conversion factor."""
    return _f(M, g) ** (1 / (g - 1))


def a_static2total(M, g):
    """Static to total sound speed conversion factor."""
    return _f(M, g) ** 0.5


def nprc(g):
    """Critical nozzle pressure ratio."""
    return p_static2total(M=1, g=g)


def compressor_s(T0_in, pi, eta_s, g_guess, table, n_iter=4):
    """Get compressor outlet flow conditions iteratively, by using the eta_s equation."""
    g = g_guess
    cp = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in * (1 + 1 / eta_s * (pi ** ((g - 1) / g) - 1))
        cp = cst.lerp_cp((T0_in + T0_out) / 2, 0)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        g = cp / (cp - cst.R_air)  # iteration update

    return (cp, T0_out, g)


def compressor_p(T0_in, pi, eta_p, g_guess, table, n_iter=4):
    """Get compressor outlet flow conditions iteratively, by using the eta_p equation."""
    g = g_guess
    cp = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in * pi ** ((g - 1) / (eta_p * g))
        cp = cst.lerp_cp((T0_in + T0_out) / 2, 0)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        g = cp / (cp - cst.R_air)  # iteration update

    return (cp, T0_out, g)


def combustion_chamber(T0_in, T0_out, eta_cc, lhv, mdot_p, far_guess, table, n_iter=4):
    """Get the fuel mass flow rate from the combustion equation, iteratively."""
    cpr_in = cst.lerp_cp((T0_in + cst.T_ref) / 2, 0)
    far = far_guess
    cpr_out = 0  # TBD
    mdot_f = far * mdot_p  # resulting mdot_f guess
    for iter in range(n_iter):
        cpr_out = cst.lerp_cp((T0_out + cst.T_ref) / 2, far)
        table.add_row(iter, cpr_out, mdot_f, far)  # keep track of iterations
        mdot_f = (
            (mdot_f + mdot_p) * cpr_out * (T0_out - cst.T_ref)
            - mdot_p * cpr_in * (T0_in - cst.T_ref)
        ) / (eta_cc * lhv)
        far = mdot_f / mdot_p

    return (mdot_f, far)


def turbine(T0_in, power, mdot_p, far, cp_guess, table, n_iter=4):
    cp = cp_guess
    g = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in - power / (mdot_p * (1 + far) * cp)
        g = cp / (cp - cst.R_air)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        cp = cst.lerp_cp((T0_in + T0_out) / 2, far)  # iteration update
    return (cp, T0_out, g)


def afterburner(T0_in, eta_ab, lhv, mdot_p, mdot_f, mdot_f_ab, table, n_iter=4):
    """Get afterburner outlet flow conditions iteratively, by using the combustion equation."""
    far = mdot_f / mdot_p
    far_ab = (mdot_f + mdot_f_ab) / mdot_p
    T0_r = cst.T_ref
    cpr_in = cst.lerp_cp((T0_in + T0_r) / 2, far)
    cpr_out = cst.lerp_cp((T0_in + T0_r) / 2, far_ab)  # initial guess
    T0_out = gr_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_r + (eta_ab * mdot_f_ab * lhv + (mdot_p + mdot_f) * cpr_in * (T0_in - T0_r)) / (
            (mdot_p + mdot_f + mdot_f_ab) * cpr_out
        )
        gr_out = cpr_out / (cpr_out - cst.R_air)
        table.add_row(iter, cpr_out, T0_out, gr_out)  # keep track of iterations
        cpr_out = cst.lerp_cp((T0_out + T0_r) / 2, far_ab)

    return cpr_out, T0_out, gr_out


def mdot_chocking(A, p0, T0, g, R=cst.R_air):
    """Choking mass flow rate."""
    return A * p0 * (g / (R * T0)) ** 0.5 * (2 / (g + 1)) ** ((g + 1) / (2 * g - 2))


def mdot_corrected(mdot, p0, T0):
    """Corrected mass flow rate."""
    return mdot * (cst.p_ref / p0) * (T0 / cst.T_ref) ** 0.5
