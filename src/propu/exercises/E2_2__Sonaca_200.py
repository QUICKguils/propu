from propu import isatmosphere as isa
from propu.constant import uconv


def main():
    print("Solving exo 2.2: Sonaca 200")

    # Statement data
    n_blades = 3
    diameter = 1750 * uconv("mm", "m")
    rho, p, T, a = isa.get_state(0)
    thrust = 1700  # [N]
    v_inf = 150 * uconv("km/hr", "m/s")

    print(f"{v_inf=} m/s")
    print(f"{p=} Pa")
