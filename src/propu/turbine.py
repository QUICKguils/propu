"""Gas turbine and jet engine analyses."""

from propu import constant as c


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


def compressor_tau(pi, eta_s, g):  # I was fed up to write this one
    """Total temperature ratio T0_out/T0_in across a compressor."""
    return 1 + 1 / eta_s * (pi ** ((g - 1) / g) - 1)


def mdot_chocking(A, p0, T0, g, R=c.R_air):
    """Choking mass flow rate."""
    return A * p0 * (g / (R * T0)) ** 0.5 * (2 / (g + 1)) ** ((g + 1) / (2 * g - 2))


def mdot_corrected(mdot, p0, T0):
    """Corrected mass flow rate."""
    return mdot * (c.p_ref / p0) * (T0 / c.T_ref) ** 0.5
