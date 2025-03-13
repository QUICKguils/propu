"""Blade element momentum techniques."""

import abc
from typing import NamedTuple

import numpy as np
from scipy import interpolate


class LocalGeometry(NamedTuple):
    span: float  # Span of one blade [m]
    n_blades: int  # Number of blades [-]
    r: float  # Radial location of the blade element centroid [m]
    dr: float  # Width of the stream tube [m]
    chord: float  # Chord [m]
    pitch: float  # Pitch [m]


class PropellerGeometry(NamedTuple):
    """Geometric description of a propeller."""

    span: float  # Span of one blade [m]
    n_blades: int  # Number of blades [-]
    stations: np.ndarray[float]  # Station locations [m]
    chords: np.ndarray[float]  # Chords evaluated at stations [m]
    pitches: np.ndarray[float]  # Pitches evaluated at stations [m]


class AirfoilPolar(NamedTuple):
    """Tabulated airfoil polar data.

    Contain the tabulated values of cl and cd, for the corresponding sampled
    angles of attack and Reynolds numbers.
    """

    aoa: np.ndarray[float]  # Sample of angles of attack [rad]
    Re: np.ndarray[float]  # Sample of Reynolds numbers [-]
    cl: np.ndarray[float]  # Corresponding lift coefficients [-]
    cd: np.ndarray[float]  # Corresponding drag coefficients [-]


class Propeller(abc.ABC):
    """Abstract base class to model a generic propeller.

    The user can create concrete propellers that derive from this base class,
    as long as it satisfies this interface. That is, as long as it has the
    following properties:
        name : str
            The name of the concrete propeller, that can be used by the creator
            to instantiate the propeller.
        airfoil: AirfoilPolar
            Tabulated airfoil polar data.
        geometry: PropellerGeometry
            Geometric description of a propeller.

    Notes
    -----
    Using a factory pattern may ba a bit overkill here, given the limited size
    of the code. However, it allows the bemt module to be uncoupled from the
    concrete propellers that can be abitrary defined elsewhere, thus making the
    module more flexible. See for example the propu.project.statement module
    for a concrete usage of the factory pattern, applied to APC propellers.

    More concretely, the author choose to use a factory pattern for the
    propellers mainly because it constitutes a good object-oriented programming
    excercise.
    """

    @property
    @abc.abstractmethod
    def name() -> str:
        pass

    @property
    @abc.abstractmethod
    def airfoil() -> AirfoilPolar:
        pass

    @property
    @abc.abstractmethod
    def geometry() -> PropellerGeometry:
        pass


class OperatingConditions(NamedTuple):
    Om: float  # Rotation speed [rad/s]
    v_inf: float  # Wind speed [m/s]
    rho: float  # Density of the air [kg/m³]
    mu: float  # Dynamic viscosity of the air [Pa*s]


class LocalBemSolution(NamedTuple):
    thrust: np.ndarray[float]  # Thrust generated in the local stream tube [N]
    torque: np.ndarray[float]  # Torque generated in the local stream tube [N*m]
    va2: np.ndarray[float]  # Local axial speed [m/s]
    vu2: np.ndarray[float]  # Local absolute tangential speed [m/s]
    wu2: np.ndarray[float]  # Local relative tangential speed [m/s]


class BemSolution(NamedTuple):
    thrust: np.ndarray[float]  # Total thrust provided by the propeller [N]
    torque: np.ndarray[float]  # Total torque provided by the propeller [N*m]
    r_dist: np.ndarray[float]  # Radial location of each blade element centroid [m]
    dr_dist: np.ndarray[float]  # Width of each blade element [m]
    thrust_dist: np.ndarray[float]  # Spanwise distribution of the stream tubes thrust [N]
    torque_dist: np.ndarray[float]  # Spanwise distribution of the stream tubes torque [N*m]
    va2_dist: np.ndarray[float]  # Spanwise distribution of va2 (= wa2) [m/s]
    vu2_dist: np.ndarray[float]  # Spanwise distribution of vu2 [m/s]
    wu2_dist: np.ndarray[float]  # Spanwise distribution of wu2 [m/s]


# FIX: crashes for v_inf = 0 m/s
def local_bem(
    lgeom: LocalGeometry, airfoil: AirfoilPolar, oper: OperatingConditions,
    *, rtol=1e-4, atol=1e-3, max_iter=50,
):
    """Momentum balance in one stream tube, for one blade.

    This function applies the blade element method to the propeller local
    geometry `lgeom`, for the specified airfoil polar data `airfoil` and
    operating conditions `oper`.
    """
    # The mass flow in null, hence va3_new and vu2p_new are NaN.

    # Initial speed guesstimates
    va3 = oper.v_inf
    vu2p = 0  # [m/s]

    # Initialize computed quantities
    thrust = 0  # [N]
    torque = 0  # [N*m]
    va2 = vu2 = wu2 = 0  # [m/s]

    # Linear interpolator methods for lift and drag coefficients
    cl_lerp = interpolate.RegularGridInterpolator((airfoil.aoa, airfoil.Re), airfoil.cl, bounds_error=False, fill_value=None)
    cd_lerp = interpolate.RegularGridInterpolator((airfoil.aoa, airfoil.Re), airfoil.cd, bounds_error=False, fill_value=None)

    # Determine forces and velocities iteratively
    has_converged = False
    for n_iter in range(1, max_iter + 1):
        # Velocity components at the propeller disk
        va2 = (oper.v_inf + va3) / 2
        vu2 = vu2p / 2
        wu2 = vu2 - oper.Om * lgeom.r
        w2 = np.sqrt(va2**2 + wu2**2)
        beta_2 = np.atan2(wu2, va2)

        aoa = lgeom.pitch - (np.pi / 2 + beta_2)
        Re = oper.rho * w2 * lgeom.chord / oper.mu

        cl = cl_lerp((aoa, Re))
        cd = cd_lerp((aoa, Re))

        lift = 0.5 * oper.rho * w2**2 * lgeom.chord * lgeom.dr * cl
        drag = 0.5 * oper.rho * w2**2 * lgeom.chord * lgeom.dr * cd

        thrust = -lgeom.n_blades * (lift * np.sin(beta_2) + drag * np.cos(beta_2))
        torque = (lgeom.n_blades * (lift * np.cos(beta_2) - drag * np.sin(beta_2)) * lgeom.r)

        # New approximations for the absolute velocity components
        mass_flow = 2 * np.pi * lgeom.r * lgeom.dr * oper.rho * va2
        va3_new = oper.v_inf + thrust / mass_flow
        vu2p_new = torque / (mass_flow * lgeom.r)

        # Stop the iterations if the speeds have sufficiently converged
        va3_ok = np.allclose(va3, va3_new, atol=atol, rtol=rtol)
        vu2p_ok = np.allclose(vu2p, vu2p_new, atol=atol, rtol=rtol)
        if va3_ok and vu2p_ok:
            has_converged = True
            break

        # Otherwise update the speeds for the next iteration
        vu2p = vu2p_new
        va3 = va3_new

    if not has_converged:
        print(f"bemt -- No convergence after {n_iter} iterations.")

    return LocalBemSolution(thrust=thrust, torque=torque, va2=va2, vu2=vu2, wu2=wu2)


def bem(prop: Propeller, oper: OperatingConditions, sdiv=20):
    """Calculate quantities derived from momentum balance of a propeller.

    This function applies the blade element method to the propeller `prop`, for
    the specified operating conditions `oper`, and the number of subdivisions
    `sdiv` in the blades.
    """
    # Virtual "cuts" in the blades, delimiting the blade elements.
    # There's no requirement to have them necessarily linearly spaced.
    r_min = prop.geometry.stations[0]
    r_max = prop.geometry.span
    rdiv = np.linspace(r_min, r_max, sdiv)

    # Radial location of the centroid of each blade element [m]
    r_dist = (rdiv[:-1] + rdiv[1:]) / 2

    # Width of each blade element [m]
    dr_dist = rdiv[1:] - rdiv[:-1]

    lerp_chords = interpolate.make_interp_spline(
        prop.geometry.stations, prop.geometry.chords, k=1
    )
    chords = lerp_chords(r_dist)

    lerp_pitches = interpolate.make_interp_spline(
        prop.geometry.stations, prop.geometry.pitches, k=1
    )
    pitches = lerp_pitches(r_dist)

    # NOTE: np.interp can not extrapolate, so the results are bad near the tip.
    # chords = np.interp(r_dist, prop.geometry.stations, prop.geometry.chords)
    # geopitches = np.interp(r_dist, prop.geometry.stations, prop.geometry.geopitches)

    thrust_dist = np.zeros(r_dist.shape)  # [N]
    torque_dist = np.zeros(r_dist.shape)  # [N*m]
    va2_dist    = np.zeros(r_dist.shape)  # [m/s]
    vu2_dist    = np.zeros(r_dist.shape)  # [m/s]
    wu2_dist    = np.zeros(r_dist.shape)  # [m/s]

    # Compute the thrust and torque for all of these locations
    for i, _ in enumerate(r_dist):
        lgeom = LocalGeometry(
            span=prop.geometry.span,
            n_blades=prop.geometry.n_blades,
            r=r_dist[i],
            dr=dr_dist[i],
            chord=chords[i],
            pitch=pitches[i],
        )

        local_bem_sol = local_bem(lgeom, prop.airfoil, oper)

        thrust_dist[i] = local_bem_sol.thrust
        torque_dist[i] = local_bem_sol.torque
        va2_dist[i]    = local_bem_sol.va2
        vu2_dist[i]    = local_bem_sol.vu2
        wu2_dist[i]    = local_bem_sol.wu2

    # Total thrust and torque exerted by the propeller
    thrust = np.sum(thrust_dist)
    torque = np.sum(torque_dist)

    return BemSolution(
        thrust=thrust,
        torque=torque,
        r_dist=r_dist,
        dr_dist=dr_dist,
        thrust_dist=thrust_dist,
        torque_dist=torque_dist,
        va2_dist=va2_dist,
        vu2_dist=vu2_dist,
        wu2_dist=wu2_dist,
    )
