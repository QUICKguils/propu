import math

from propu import constant as c


def main():
    print("Solving exo 3.1: Comparison of compressors")

    # Statement data

    class A:
        p0_ratio = 4.5  # Total pressure ratio
        eta_s = 0.85  # Isentropic efficiency

    class B:
        p_2 = 90 * c.uconv("psi", "Pa")  # Outlet pressure
        T0_2 = c.mconv(480, "degF", "K")  # Outlet total temperature
        v_2 = 500 * c.uconv("ft/s", "m/s")  # Outlet velocity

    p_1 = 14.7 * c.uconv("psi", "Pa")  # Inlet pressure
    T_1 = c.mconv(60, "degF", "K")  # Inlet temperature
    g = c.gamma_air  # Adiabatic index

    # Resolution
    # NOTE: see complement (3)

    # Compressor A.
    # Just revert the formula (3.32), p. 26-2
    A.eta_p = ((g - 1) * math.log(A.p0_ratio)) / (
        g * math.log((A.p0_ratio ** ((g - 1) / g) - 1) / (A.eta_s) + 1)
    )

    # Compressor B.
    # Calculate T2, from its def. (B.8),
    # then find eta_p from (3.29)
    B.T2 = B.T0_2 - (B.v_2**2 * (g - 1)) / (2 * c.R_air * g)
    B.eta_p = (g - 1) / g * math.log(B.p_2 / p_1) / math.log(B.T2 / T_1)

    print(f"Polytropic efficiency for compressor A: {A.eta_p:.2%}")
    print(f"Polytropic efficiency for compressor B: {B.eta_p:.2%}")
