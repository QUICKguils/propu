"""isatmosphere -- International Standard Atmosphere.

Reference
---------
[U.S. Standard Atmosphere, 1976, U.S. Government Printing Office, Washington, D.C., 1976]
(https://en.wikipedia.org/wiki/U.S._Standard_Atmosphere)

| Layer | ``h`` (m) | ``p`` (Pa) | ``T`` (K) | ``α`` (K/m) |
| ----- | -----     | -------    | ------    | ----------- |
| 0     | 0         | 101325     | 288.15    | -0.0065     |
| 1     | 11000     | 22632.1    | 216.65    | 0           |
| 2     | 20000     | 5474.89    | 216.65    | 0.001       |
| 3     | 32000     | 868.019    | 228.65    | 0.0028      |
| 4     | 47000     | 110.906    | 270.65    | 0           |
| 5     | 51000     | 66.9389    | 270.65    | -0.0028     |
| 6     | 71000     | 3.95642    | 214.65    | -0.002      |

Credits
-------
The main `get_state` function is mainly a python rewrite of:
https://github.com/AlexS12/FlightMechanicsUtils.jl/blob/main/src/atmosphere.jl
"""

import numpy as np

from ..util import constant as c
from ..util.constant import ur

def get_state(height: ur.Quantity):
    """(rho, p, T, a) = get_state(height)

    Get the density, pressure, temperature and sound speed for the given
    geopotential height, according to the International Standard Atmosphere.
    """
    if 0*ur.m <= height < 11_000*ur.m:  # Troposphere
        alpha  = -0.0065 * ur("K/m")
        T0 = 288.15 * ur.K
        p0 = 101_325.0 * ur.Pa

        T = T0 + alpha * height
        p = p0 * (T0 / T)**(1*ur.g_n / (c.R_air * alpha))

    elif 11_000*ur.m <= height < 20_000*ur.m:  # Tropopause
        T  = 216.65 * ur.K
        p0 = 22_632.1 * ur.Pa
        h0 = 11_000*ur.m

        p = p0 * np.exp(-1*ur.g_n * (height - h0) / (c.R_air * T))

    elif 20_000*ur.m <= height < 32_000*ur.m:  # Stratosphere 1
        alpha  = 0.001 * ur("K/m")
        T0 = 216.65 * ur.K
        p0 = 5474.89 * ur.Pa
        h0 = 20_000*ur.m

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T)**(1*ur.g_n / (c.R_air * alpha))

    elif 32_000*ur.m <= height < 47_000*ur.m:  # Stratosphere 2
        alpha  = 0.0028 * ur("K/m")
        T0 = 228.65 * ur.K
        p0 = 868.019 * ur.Pa
        h0 = 32_000*ur.m

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T)**(1*ur.g_n / (c.R_air * alpha))

    elif 47_000*ur.m <= height < 51_000*ur.m:  # Stratopause
        T  = 270.65 * ur.K
        p0 = 110.906 * ur.Pa
        h0 = 47_000*ur.m

        p = p0 * np.exp(-1*ur.g_n * (height - h0) / (c.R_air * T))

    elif 51_000*ur.m <= height < 71_000*ur.m:  # Mesosphere 1
        alpha  = -0.0028 * ur("K/m")
        T0 = 270.65 * ur.K
        p0 = 66.9389 * ur.Pa
        h0 = 51_000*ur.m

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T)**(1*ur.g_n / (c.R_air * alpha))

    elif 71_000*ur.m <= height <= 84_500*ur.m:  # Mesosphere 2
        alpha  = -0.002 * ur("K/m")
        T0 = 214.65 * ur.K
        p0 = 3.95642 * ur.Pa
        h0 = 71_000 * ur.m

        T = T0 + alpha * (height - h0)
        p = p0 * (T0 / T)**(1*ur.g_n / (c.R_air * alpha))

    else:
        raise AttributeError("height must be between 0 m and 84500 m")

    rho = p / (c.R_air * T)
    a = np.sqrt(c.gamma_air * c.R_air * T)

    return (
        rho.to("kg/m**3"),
        p.to(ur.Pa),
        T.to(ur.K),
        a.to("m/s")
    )
