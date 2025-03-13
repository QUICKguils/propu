"""Define the general project statement data."""

import pathlib

import numpy as np

from propu import bemt
from propu.constant import uconv

_PROJECT_PATH = pathlib.Path(__file__).parent
_DATA_PATH = _PROJECT_PATH / "res"

airfoil_profile = np.loadtxt(str(_DATA_PATH / "clarky.dat"), skiprows=1)
airfoil_polar = bemt.AirfoilPolar(
    # XXX: see clarkypolarsRe.py: maybe transform aoa with arctan2
    aoa = np.loadtxt(str(_DATA_PATH / "aoa.dat")),  # NOTE: aoa is given in [rad]
    Re  = np.loadtxt(str(_DATA_PATH / "reynolds.dat")),
    cl  = np.loadtxt(str(_DATA_PATH / "cl.dat")),
    cd  = np.loadtxt(str(_DATA_PATH / "cd.dat")),
)


class APC_9x45(bemt.Propeller):
    name = "apc_9x45"

    airfoil = airfoil_polar

    def _geometry():
        span = (9 / 2) * uconv("in", "m")
        n_blades = 2
        datafile = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x45_geom.txt"), skiprows=1)
        stations = datafile[:, 0] * span
        chords = datafile[:, 1] * span
        pitches = np.deg2rad(datafile[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_9x6(bemt.Propeller):
    name = "apc_9x6"

    airfoil = airfoil_polar

    def _geometry():
        span = (9 / 2) * uconv("in", "m")
        n_blades = 2
        datafile = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x6_geom.txt"), skiprows=1)
        stations = datafile[:, 0] * span
        chords = datafile[:, 1] * span
        pitches = np.deg2rad(datafile[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_11x7(bemt.Propeller):
    name = "apc_11x7"

    airfoil = airfoil_polar

    def _geometry():
        span = (11 / 2) * uconv("in", "m")
        n_blades = 2
        datafile = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x7_geom.txt"), skiprows=1)
        stations = datafile[:, 0] * span
        chords = datafile[:, 1] * span
        pitches = np.deg2rad(datafile[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_11x10(bemt.Propeller):
    name = "apc_11x10"

    airfoil = airfoil_polar

    def _geometry():
        span = (11 / 2) * uconv("in", "m")
        n_blades = 2
        datafile = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x10_geom.txt"), skiprows=1)
        stations = datafile[:, 0] * span
        chords = datafile[:, 1] * span
        pitches = np.deg2rad(datafile[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


APC_dict = {
    APC_9x45.name: APC_9x45,
    APC_9x6.name: APC_9x6,
    APC_11x7.name: APC_11x7,
    APC_11x10.name: APC_11x10,
}
"""Holds refences to all the concrete propellers defined in this module."""


def create_apc_propeller(propeller_name: str) -> bemt.Propeller:
    """Concrete propeller creator."""
    return APC_dict[propeller_name]()
