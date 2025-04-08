"""Define the general project statement data."""

import pathlib

import numpy as np

from propu import bemt
from propu.constant import uconv

_PROJECT_PATH = pathlib.Path(__file__).parent
_DATA_PATH = _PROJECT_PATH / "res"

airfoil = bemt.Airfoil(
    profile = np.loadtxt(str(_DATA_PATH / "clarky.dat"), skiprows=1),
    aoa = np.loadtxt(str(_DATA_PATH / "aoa.dat")),
    Re  = np.loadtxt(str(_DATA_PATH / "reynolds.dat")),
    cl  = np.loadtxt(str(_DATA_PATH / "cl.dat")),
    cd  = np.loadtxt(str(_DATA_PATH / "cd.dat")),
)


class APC_9x45(bemt.Propeller):
    keyword = "apce_9x45"
    pretty_name = "APC 9X4.5"

    airfoil = airfoil

    def _geometry():
        span = (9 / 2) * uconv("in", "m")
        n_blades = 2
        data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x45_geom.txt"), skiprows=1)
        stations = data[:, 0] * span
        chords = data[:, 1] * span
        pitches = np.deg2rad(data[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_9x6(bemt.Propeller):
    keyword = "apce_9x6"
    pretty_name = "APC 9X6"

    airfoil = airfoil

    def _geometry():
        span = (9 / 2) * uconv("in", "m")
        n_blades = 2
        data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x6_geom.txt"), skiprows=1)
        stations = data[:, 0] * span
        chords = data[:, 1] * span
        pitches = np.deg2rad(data[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_11x7(bemt.Propeller):
    keyword = "apce_11x7"
    pretty_name = "APC 11X7"

    airfoil = airfoil

    def _geometry():
        span = (11 / 2) * uconv("in", "m")
        n_blades = 2
        data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x7_geom.txt"), skiprows=1)
        stations = data[:, 0] * span
        chords = data[:, 1] * span
        pitches = np.deg2rad(data[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


class APC_11x10(bemt.Propeller):
    keyword = "apce_11x10"
    pretty_name = "APC 11X10"

    airfoil = airfoil

    def _geometry():
        span = (11 / 2) * uconv("in", "m")
        n_blades = 2
        data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x10_geom.txt"), skiprows=1)
        stations = data[:, 0] * span
        chords = data[:, 1] * span
        pitches = np.deg2rad(data[:, 2])
        return bemt.PropellerGeometry(
            span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches,
        )

    geometry = _geometry()


APC_dict = {
    APC_9x45.keyword: APC_9x45,
    APC_9x6.keyword: APC_9x6,
    APC_11x7.keyword: APC_11x7,
    APC_11x10.keyword: APC_11x10,
}
"""Holds refences to all the concrete propellers defined in this module."""


def create_apc_propeller(prop_keyword: str) -> bemt.Propeller:
    """Concrete propeller creator."""
    return APC_dict[prop_keyword]()
