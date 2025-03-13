"""International Standard Atmosphere."""

import numpy as np

from propu import constant as c
from propu.constant import uconv, ureg


def get_state(height: float) -> tuple[float, float, float, float]:
    """Get the atmospheric state variables at the given height.

    Get the density, pressure, temperature and sound speed for the given
    geopotential height, according to the International Standard Atmosphere [1].

    Parameters
    ----------
    height : float
        The geopotential height (m).
        Must be between 0 m and 84_500 m.

    Returns
    -------
    rho : float
        Density (kg/m³).
    p : float
        Pressure (Pa).
    T : float
        Temperature (K).
    a : float
        Sound speed (m/s).

    Raises
    ------
    ValueError
        If `height` is not between 0 m and 84_500 m.

    Notes
    -----
    This function is mainly a python rewrite of:
    https://github.com/AlexS12/FlightMechanicsUtils.jl/blob/main/src/atmosphere.jl

    References
    ----------
    .. [1] U.S. Standard Atmosphere,
       1976, U.S. Government Printing Office, Washington, D.C., 1976
       (https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere)

    ===== ======= ======== ======= =========
    Layer `h` (m) `p` (Pa) `T` (K) `α` (K/m)
    ===== ======= ======== ======= =========
        0       0  101325   288.15   -0.0065
        1   11000  22632.1  216.65   0
        2   20000  5474.89  216.65   0.001
        3   32000  868.019  228.65   0.0028
        4   47000  110.906  270.65   0
        5   51000  66.9389  270.65   -0.0028
        6   71000  3.95642  214.65   -0.002
    ===== ======= ======== ======= =========
    """
    if 0 <= height < 11_000:  # Troposphere
        alpha = -0.0065  # [K/m]
        T0 = 288.15  # [K]
        p0 = 101_325.0  # [Pa]

        T = T0 + alpha * height
        p = p0 * (T0 / T) ** (uconv("g_n", "m/s**2") / (c.R_air * alpha))

    elif 11_000 <= height < 20_000:  # Tropopause
        T = 216.65  # [K]
        p0 = 22_632.1  # [Pa]
        h0 = 11_000  # [m]

        p = p0 * np.exp(-1 * ureg.g_n * (height - h0) / (c.R_air * T))

    elif 20_000 <= height < 32_000:  # Stratosphere 1
        alpha = 0.001  # [K/m]
        T0 = 216.65  # [K]
        p0 = 5474.89  # [Pa]
        h0 = 20_000  # [m]

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T) ** (uconv("g_n", "m/s**2") / (c.R_air * alpha))

    elif 32_000 <= height < 47_000:  # Stratosphere 2
        alpha = 0.0028  # [K/m]
        T0 = 228.65  # [K]
        p0 = 868.019  # [Pa]
        h0 = 32_000  # [m]

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T) ** (uconv("g_n", "m/s**2") / (c.R_air * alpha))

    elif 47_000 <= height < 51_000:  # Stratopause
        T = 270.65  # [K]
        p0 = 110.906  # [Pa]
        h0 = 47_000  # [m]

        p = p0 * np.exp(-uconv("g_n", "m/s**2") * (height - h0) / (c.R_air * T))

    elif 51_000 <= height < 71_000:  # Mesosphere 1
        alpha = -0.0028  # [K/m]
        T0 = 270.65  # [K]
        p0 = 66.9389  # [Pa]
        h0 = 51_000  # [m]

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T) ** (uconv("g_n", "m/s**2") / (c.R_air * alpha))

    elif 71_000 <= height <= 84_500:  # Mesosphere 2
        alpha = -0.002  # [K/m]
        T0 = 214.65  # [K]
        p0 = 3.95642  # [Pa]
        h0 = 71_000  # [m]

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T) ** (uconv("g_n", "m/s**2") / (c.R_air * alpha))

    else:
        raise ValueError("height must be between 0m and 84_500m")

    rho = p / (c.R_air * T)
    a = np.sqrt(c.gamma_air * c.R_air * T)

    state = (rho, p, T, a)
    return tuple(float(s) for s in state)
