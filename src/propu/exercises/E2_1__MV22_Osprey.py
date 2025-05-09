import numpy as np

from propu.constant import uconv
from propu.isatmosphere import get_state


def main():
    print("Solving exo 2.1: MV-22 Osprey")

    # Statement data
    disk_loading = 11.4 * uconv("psf", "Pa")
    v_inf = 80 * uconv("knot", "m/s")
    rho_sl = get_state(0).rho

    # Resolution
    # NOTE: this commented eta_p is not exact (see personnal notes)
    # eta_p = 2v_inf / (2v_inf + (disk_loading/(rho_sl*v_inf)))
    eta_p = 1 / (0.5 + np.sqrt(0.25 + disk_loading / (2 * rho_sl * v_inf**2)))

    print(f"{eta_p=}")
