"""Jet engine cycle analyses."""

from typing import NamedTuple
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


def get_nprc(g):
    """Critical nozzle pressure ratio."""
    return p_static2total(M=1, g=g)


def compressor_s(T0_in, pi, eta_s, g_guess, table: IterTable, n_iter=4):
    """Get compressor outlet flow conditions iteratively, by using the eta_s equation."""
    g = g_guess
    cp = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in * (1 + 1 / eta_s * (pi ** ((g - 1) / g) - 1))
        cp = cst.lerp_cp((T0_in + T0_out) / 2, 0)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        g = cp / (cp - cst.R_air)  # iteration update

    return (cp, T0_out, g)


def compressor_p(T0_in, pi, eta_p, g_guess, table: IterTable, n_iter=4):
    """Get compressor outlet flow conditions iteratively, by using the eta_p equation."""
    g = g_guess
    cp = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in * pi ** ((g - 1) / (eta_p * g))
        cp = cst.lerp_cp((T0_in + T0_out) / 2, 0)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        g = cp / (cp - cst.R_air)  # iteration update

    return (cp, T0_out, g)


def combustion_chamber(T0_in, T0_out, eta_cc, lhv, mdot_p, far_guess, table: IterTable, n_iter=4):
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


def turbine(T0_in, power, mdot_p, far, cp_guess, table: IterTable, n_iter=4):
    """Get turbine outlet flow conditions iteratively, by using the turbine power equation."""
    cp = cp_guess
    g = T0_out = 0  # TBD
    for iter in range(n_iter):
        T0_out = T0_in - power / (mdot_p * (1 + far) * cp)
        g = cp / (cp - cst.R_air)
        table.add_row(iter, cp, T0_out, g)  # keep track of iterations
        cp = cst.lerp_cp((T0_in + T0_out) / 2, far)  # iteration update
    return (cp, T0_out, g)


def afterburner(T0_in, eta_ab, lhv, mdot_p, mdot_f, mdot_f_ab, table: IterTable, n_iter=4):
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

    return (cpr_out, T0_out, gr_out)


def nozzle_sonic(T0, p0, far, table: IterTable, n_iter=4):
    """Get nozzle sonic flow conditions iteratively, by using isentropic relations."""
    M = 1  # By def. of sonic conditions

    cp = cst.lerp_cp(T0, far)  # initial guess
    nprc = T = g = 0  # TBD
    for iter in range(n_iter):
        g = cp / (cp - cst.R_air)
        nprc = get_nprc(g)  # critical nozzle pressure ratio
        T = T0 / T_static2total(M=M, g=g)  # Expansion/compression to sonic conditions
        table.add_row(iter, cp, T, g, nprc)  # keep track of iterations
        cp = cst.lerp_cp((T + T0) / 2, far)  # iteration update

    p = p0 / nprc
    rho = p / (cst.R_air * T)
    a = (g * cst.R_air * T) ** 0.5  # = speed, as sonic conditions (M = 1)

    class SonicConditions(NamedTuple):
        rho: float  # density
        p: float  # pressure
        T: float  # temperature
        a: float  # sonic speed
        cp: float  # heat capacity at constant pressure
        g: float  # adiabatic index
        nprc: float  # critical nozzle pressure ratio

    return SonicConditions(rho=rho, p=p, T=T, a=a, cp=cp, g=g, nprc=nprc)


def nozzle_adapted(T0, npr, far, g_guess, table: IterTable, n_iter=4):
    """Get adapted nozzle exhaust flow conditions iteratively, by using isentropic relations."""
    g = g_guess  # initial guess
    M = T = 0  # TBD
    for iter in range(n_iter):
        T = T0 * npr ** ((1 - g) / g)  # Isentropic relations
        cp = cst.lerp_cp((T0 + T) / 2, far)
        M = (((T0 / T) - 1) * 2 / (g - 1)) ** 0.5
        table.add_row(iter, cp, T, g, M)  # keep track of iterations
        g = cp / (cp - cst.R_air)  # iteration update

    a = (g * cst.R_air * T) ** 0.5
    v = M * a

    class AdaptedExhaustConditions(NamedTuple):
        T: float  # temperature
        a: float  # sonic speed
        M: float  # Mach number
        v: float  # speed
        cp: float  # heat capacity at constant pressure
        g: float  # adiabatic index

    return AdaptedExhaustConditions(T=T, a=a, M=M, v=v, cp=cp, g=g)


def mdot_chocking(A, p0, T0, g, R=cst.R_air):
    """Choking mass flow rate."""
    return A * p0 * (g / (R * T0)) ** 0.5 * (2 / (g + 1)) ** ((g + 1) / (2 * g - 2))


def mdot_corrected(mdot, p0, T0):
    """Corrected mass flow rate."""
    return mdot * (cst.p_ref / p0) * (T0 / cst.T_ref) ** 0.5
