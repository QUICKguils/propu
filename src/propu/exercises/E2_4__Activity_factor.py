import numpy as np

from propu.constant import uconv


def main():
    print("Solving exo 2.4: Activity Factor")

    # Statement data
    D = 8 * uconv("in", "m")
    stations = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]) * D/2
    apc_8x6S = np.array([13.5, 12.6, 17.6, 21.2, 23.5, 24.4, 23.3, 20.8, 16.0, 6.3]) / 1e3
    apc_8x6 = np.array([19.2, 19.0, 18.8, 18.8, 18.6, 18.6, 18.2, 14.8, 12.5, 7.1]) / 1e3
    apc_8x6E = np.array([16.5, 13.1, 19.6, 24.2, 24.6, 23.3, 19.4, 14.4, 11.5, 2.4]) / 1e3

    # Resolution
    # Just apply def. of activity factor, p.10-2
    def af(c_dist):
        return 1e5 / D**5 * np.trapezoid(c_dist * stations**3, stations)

    print(f"AF for apc8x6s: {af(apc_8x6S):.2f}")
    print(f"AF for apc8x6: {af(apc_8x6):.2f}")
    print(f"AF for apc8x6e: {af(apc_8x6E):.2f}")
