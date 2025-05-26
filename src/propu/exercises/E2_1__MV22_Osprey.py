import propu.constant as cst


def main():
    print("Solving exo 2.1: MV-22 Osprey")
    # NOTE: see complement (1)

    # Statement data
    dl = cst.mconv(11.4, "psf", "Pa")  # Disk loading
    v_inf = cst.mconv(80, "knot", "m/s")
    rho = cst.get_isa(0).rho

    # Resolution
    v_e = (2*dl/rho + v_inf**2) ** 0.5
    eta_p = 2*v_inf/(v_inf+v_e)

    print(f"Propulsive efficiency: {100*eta_p:2.2f} %")
