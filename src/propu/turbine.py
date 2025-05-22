"""Gas turbine analyses."""

from propu import constant as c


def T_static2total(M, gamma=c.gamma_air):
    """Static to total temperature conversion factor."""
    return 1 + (gamma - 1) / 2 * M**2


def p_static2total(M, gamma=c.gamma_air):
    """Static to total pressure conversion factor."""
    return T_static2total(M, gamma) ** (gamma / (gamma - 1))


def rho_static2total(M, gamma=c.gamma_air):
    """Static to total density conversion factor."""
    return T_static2total(M, gamma) ** (1 / (gamma - 1))


def a_static2total(M, gamma=c.gamma_air):
    """Static to total sound speed conversion factor."""
    return T_static2total(M, gamma) ** 0.5


def mdot_chocking(section_area, p_total, T_total, gamma=c.gamma_air, R=c.R_air):
    """Choking mass flow rate."""
    return (
        section_area
        * p_total
        * (gamma / (R * T_total)) ** 0.5
        * (2 / (gamma + 1)) ** ((gamma + 1) / (2 * gamma - 2))
    )


def mdot_corrected(pt, Tt, m_flow):
    """Corrected/reduced mass flow."""
    return m_flow * (c.p_ref / pt) * (Tt / c.T_ref) ** 0.5
