"""Define the general project statement data."""

import pathlib

import numpy as np

from propu import bemt
from propu.constant import uconv
from propu.isatmosphere import get_state

_PROJECT_PATH = pathlib.Path(__file__).parent
_DATA_PATH = _PROJECT_PATH / "res"

# Conventional air properties (STP, sea level)
rho = get_state(0).rho  # Density [kg/m³]
mu = 17.89e-6  # Dynamic viscosity [Pa*s]

airfoil = bemt.Airfoil(
    profile = np.loadtxt(str(_DATA_PATH / "clarky" / "profile.dat"), skiprows = 1),
    aoa     = np.loadtxt(str(_DATA_PATH / "clarky" / "aoa.dat")),
    Re      = np.loadtxt(str(_DATA_PATH / "clarky" / "reynolds.dat")),
    cl      = np.loadtxt(str(_DATA_PATH / "clarky" / "cl.dat")),
    cd      = np.loadtxt(str(_DATA_PATH / "clarky" / "cd.dat")),
)


def _create_APC_9x45() -> bemt.Propeller:
    span = (9 / 2) * uconv("in", "m")
    n_blades = 2
    data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x45_geom.txt"), skiprows=1)
    stations = data[:, 0] * span
    chords = data[:, 1] * span
    pitches = np.deg2rad(data[:, 2])
    geometry = bemt.PropellerGeometry(
        span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches
    )
    return bemt.Propeller(
        keyword="apce_9x45", pretty_name="APC 9x4.5", airfoil=airfoil, geometry=geometry
    )


def _create_APC_9x6() -> bemt.Propeller:
    span = (9 / 2) * uconv("in", "m")
    n_blades = 2
    data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_9x6_geom.txt"), skiprows=1)
    stations = data[:, 0] * span
    chords = data[:, 1] * span
    pitches = np.deg2rad(data[:, 2])
    geometry = bemt.PropellerGeometry(
        span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches
    )
    return bemt.Propeller(
        keyword="apce_9x6", pretty_name="APC 9x6", airfoil=airfoil, geometry=geometry
    )


def _create_APC_11x7() -> bemt.Propeller:
    span = (11 / 2) * uconv("in", "m")
    n_blades = 2
    data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x7_geom.txt"), skiprows=1)
    stations = data[:, 0] * span
    chords = data[:, 1] * span
    pitches = np.deg2rad(data[:, 2])
    geometry = bemt.PropellerGeometry(
        span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches
    )
    return bemt.Propeller(
        keyword="apce_11x7", pretty_name="APC 11x7", airfoil=airfoil, geometry=geometry
    )


def _create_APC_11x10() -> bemt.Propeller:
    span = (11 / 2) * uconv("in", "m")
    n_blades = 2
    data = np.loadtxt(str(_DATA_PATH / "apc" / "apce_11x10_geom.txt"), skiprows=1)
    stations = data[:, 0] * span
    chords = data[:, 1] * span
    pitches = np.deg2rad(data[:, 2])
    geometry = bemt.PropellerGeometry(
        span=span, n_blades=n_blades, stations=stations, chords=chords, pitches=pitches
    )
    return bemt.Propeller(
        keyword="apce_11x10", pretty_name="APC 11x10", airfoil=airfoil, geometry=geometry
    )


propellers = [
    _create_APC_9x45(),
    _create_APC_9x6(),
    _create_APC_11x7(),
    _create_APC_11x10(),
]
"""Instance of all the propellers defined in this module."""


def get_apc_propeller(keyword: str) -> bemt.Propeller:
    """Get the desired APC propeller."""
    propeller = next((prop for prop in propellers if prop.keyword == keyword), None)
    if propeller is None:
        raise ValueError("Bro pls gimme a valid propeller keyword.")
    return propeller
