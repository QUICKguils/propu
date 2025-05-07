"""Code developed to answer the fourth part of the project.

For the latter configuration, investigate strategies to maintain a constant
rotational velocity of 8000 RPM with a fixed power input of 0.31 hp and varying
advance velocity.
"""

from typing import NamedTuple

import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state
from propu.project import statement as stm


class Performance(NamedTuple):
    J: np.ndarray[float]
    eta: np.ndarray[float]
    drag: np.ndarray[float]
    converged: np.ndarray[bool]


class Solution(NamedTuple):
    """Associate a propeller to its computed performances."""
    prop: bemt.Propeller
    perf: Performance


def main(*, out_enabled=True) -> Solution:
    """Execute the fourth part of the project."""
    prop = stm.APC_11x7()  # As determined in part3
    perf = compute_performance(prop)
    sol = Solution(prop=prop, perf=perf)

    if out_enabled:
        plot_solution(sol)

    return sol

def get_operating_conditions(prop: bemt.Propeller, J: float) -> bemt.OperatingConditions:
    Om = 8000 * uconv("rpm", "rad/s")
    D = 2 * prop.geometry.span
    rho, _, _, _ = get_state(0)  # Density of the air at sea level [kg/m³]
    mu = 17.89e-6  # Dynamic viscosity of the air [Pa*s]
    v_inf = J * D * Om * uconv("rad/s", "rps")
    return bemt.OperatingConditions(Om=Om, v_inf=v_inf, rho=rho, mu=mu)


def compute_performance(prop: bemt.Propeller) -> Performance:
    POWER_IMPOSED = 0.31 * uconv("hp", "W")
    return Performance()


def plot_solution(sol: Solution) -> None:
    """Plot the computed solutions.

    The matplotlib object instances are kept as closures,
    to ensure their idempotence through multiple plot calls.
    """
    pass
