import numpy as np

from propu.constant import uconv, get_isa


def main():
    print("Solving exo 2.1: MV-22 Osprey")

    # Statement data
    l = 11.4 * uconv("psf", "Pa")  # Disk loading
    v_inf = 80 * uconv("knot", "m/s")
    rho = get_isa(0).rho

    # Resolution
    # NOTE: see complement (1)
    v_e = np.sqrt(2*l/rho + v_inf**2)
    eta_p = 2*v_inf/(v_inf+v_e)

    print(f"Propulsive efficiency: {100*eta_p:2.2f} %")
