"""Gas turbine analyses."""

import numpy as np

from propu import constant as c


def T_static2total(M, gamma):
    """Static to total temperature factor from Mach number and heat capacity ratio."""
    return 1 + (gamma - 1) / 2 * M**2


def p_static2total(M, gamma):
    """Static to total pressure factor from Mach number and heat capacity ratio."""
    return T_static2total(M, gamma) ** (gamma / (gamma - 1))


def mdot_chocking(pt, Tt, R, gamma, A):
    """Choking mass flow."""
    return pt / np.sqrt(gamma*R*Tt) * (2/(gamma+1)) ** ((gamma+1)/2/(gamma-1)) * A


def mdot_corrected(pt, Tt, m_flow):
    """Corrected/reduced mass flow."""
    return m_flow * (c.p_ref / pt) * np.sqrt(Tt / c.T_ref)
