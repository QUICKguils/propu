import numpy as np

from propu.constant import uconv
from propu.isatmosphere import get_state


def main():
    print("Solving exo 2.2: Sonaca 200")

    # Statement data
    D = 1750 * uconv("mm", "m")  # Diameter
    air = get_state(0)  # ISA properties at sea level
    T = 1700  # [N]
    vinf = 150 * uconv("km/hr", "m/s")  # Upstream wind speed

    # Resolution
    Ap = np.pi * D**2 / 4
    ve = np.sqrt(2 * T / (air.rho * Ap) + vinf**2)  # Similar to E2_1
    v2 = (vinf + ve) / 2
    eta = 2 * vinf / (vinf + ve)
    pdiff = 0.5 * air.rho * (v2**2 - vinf**2)  # Just Bernouilli

    print(f"Front propeller speed: {v2:.2f} m/s")
    print(f"Propulsive efficiency: {100 * eta:.2f} %")
    print(f"Pressure difference: {pdiff:.2f} Pa")
    print(f"Final wake speed: {ve:.2f} m/s")
