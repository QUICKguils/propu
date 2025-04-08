"""Blade element momentum techniques."""

import abc
from typing import NamedTuple

import numpy as np
from scipy import interpolate


class Airfoil(NamedTuple):
    """Tabulated Airfoil data.

    Contain the tabulated coordinates of the airfoil profile, as well as
    the tabulated values of cl and cd, for the corresponding sampled angles of
    attack and Reynolds numbers.
    """
    profile: np.ndarray[float]  # Sample of airfoil coordinates [m]
    aoa: np.ndarray[float]  # Sample of angles of attack [rad]
    Re: np.ndarray[float]  # Sample of Reynolds numbers [-]
    cl: np.ndarray[float]  # Corresponding lift coefficients [-]
    cd: np.ndarray[float]  # Corresponding drag coefficients [-]


class PropellerGeometry(NamedTuple):
    """Geometric description of a propeller."""
    span: float  # Span of one blade [m]
    n_blades: int  # Number of blades [-]
    stations: np.ndarray[float]  # Station locations [m]
    chords: np.ndarray[float]  # Chords evaluated at stations [m]
    pitches: np.ndarray[float]  # Pitches evaluated at stations [m]


class Propeller(abc.ABC):
    """Abstract base class to model a generic propeller.

    The user can create concrete propellers that derive from this base class,
    as long as it satisfies this interface. That is, as long as it has the
    following properties:
        keyword : str
            The reference name of the concrete propeller. This can be used as a
            keyword by the creator, to instantiate the propeller.
        pretty_name: str
            The well-formatted name of the concrete propeller. This can be
            used, for example, when printing results or in plot legends.
        airfoil: AirfoilPolar
            Tabulated airfoil polar data.
        geometry: PropellerGeometry
            Geometric description of a propeller.

    Notes
    -----
    Using a factory pattern may be a bit overengineered here, given the limited
    size of the code. However, it allows the bemt module to be uncoupled from
    the concrete propellers that can be abitrary defined elsewhere, thus making
    this module more flexible.
    See for example the propu.project.statement module for a concrete usage of
    this factory pattern, applied to APC propellers.

    More concretely, the author choose to use a factory pattern for the
    propellers mainly because it constitutes a good object-oriented programming
    excercise.
    """
    @property
    @abc.abstractmethod
    def keyword() -> str:
        pass

    @property
    @abc.abstractmethod
    def pretty_name() -> str:
        pass

    @property
    @abc.abstractmethod
    def airfoil() -> Airfoil:
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


class Location(NamedTuple):
    """Radial location of one blade element."""
    r: float  # Radial location of the blade element centroid [m]
    dr: float  # Width of the stream tube [m]


class LocalBemSolution(NamedTuple):
    # Context
    prop: Propeller  # Which propeller is used ?
    oper: OperatingConditions  # Under what operating conditions ?
    loc: Location  # At which location on the blade (blade element) ?
    # Results
    va2: float  # Local axial speed [m/s]
    vu2: float  # Local absolute tangential speed [m/s]
    wu2: float  # Local relative tangential speed [m/s]
    thrust: float  # Thrust generated in the local stream tube [N]
    torque: float  # Torque generated in the local stream tube [N*m]
    # Status
    flag_converged: bool  # Raised if the solution has converged.


class BemSolution:
    def __init__(self, lsols: np.ndarray[LocalBemSolution]):
        # Context
        self.prop = lsols[0].prop  # Which propeller is used ?
        self.oper = lsols[0].oper  # Under what operating conditions ?
        # Extract data from local bem solutions in separate lists
        self.r_dist = np.array([sol.loc.r for sol in lsols])
        self.dr_dist = np.array([sol.loc.dr for sol in lsols])
        self.flag_converged_dist = np.array([sol.flag_converged for sol in lsols])
        self.va2_dist = np.array([sol.va2 for sol in lsols])
        self.vu2_dist = np.array([sol.vu2 for sol in lsols])
        self.wu2_dist = np.array([sol.wu2 for sol in lsols])
        self.thrust_dist = np.array([sol.thrust for sol in lsols])
        self.torque_dist = np.array([sol.torque for sol in lsols])
        self.power_dist = self.torque_dist * self.oper.Om
        # Build global propeller performances quantities
        self.thrust = np.sum(self.thrust_dist)  # Total thrust provided by the propeller [N]
        self.torque = np.sum(self.torque_dist)  # Total torque provided by the propeller [N*m]
        self.power = np.sum(self.power_dist)  # Total power neeed to drive the propeller [W]
        self.eta = self.thrust/self.power * self.oper.v_inf  # Propulsive efficiency [-]


def local_bem(
    prop: Propeller, oper: OperatingConditions, loc: Location,
    *, rtol=1e-5, atol=1e-5, max_iter=100,
) -> LocalBemSolution:
    """Momentum balance in one stream tube, for one blade.

    This function applies the blade element method to the specified propeller
    `prop`, under the operating conditions `oper` and for the balde element
    located at `loc`.
    """
    # Initial speed guesstimates
    #
    # NOTE: Avoid to set null speed for va3.
    # If this is the case, then the computed mass flow in null,
    # hence va3_new and vu2p_new are NaN.
    va3 = max(oper.v_inf, 1)
    vu2p = 0  # [m/s]

    # Linear interpolation on the blade span to retrieve pitch and chord at `loc`
    chord_lerp = interpolate.make_interp_spline(prop.geometry.stations, prop.geometry.chords, k=1)
    pitch_lerp = interpolate.make_interp_spline(prop.geometry.stations, prop.geometry.pitches, k=1)
    chord = chord_lerp(loc.r)
    pitch = pitch_lerp(loc.r)

    # Linear interpolator methods for lift and drag coefficients
    cl_lerp = interpolate.RegularGridInterpolator(
        (prop.airfoil.aoa, prop.airfoil.Re), prop.airfoil.cl,
        bounds_error=False, fill_value=None
    )
    cd_lerp = interpolate.RegularGridInterpolator(
        (prop.airfoil.aoa, prop.airfoil.Re), prop.airfoil.cd,
        bounds_error=False, fill_value=None
    )

    # Determine forces and velocities iteratively
    relax = 0.3
    flag_converged = False
    for n_iter in range(1, max_iter + 1):
        # Velocity components at the propeller disk
        va2 = (oper.v_inf + va3) / 2
        vu2 = vu2p / 2
        wu2 = vu2 - oper.Om * loc.r
        w2 = np.sqrt(va2**2 + wu2**2)
        beta_2 = np.atan2(wu2, va2)

        aoa = pitch - (np.pi / 2 + beta_2)
        Re = oper.rho * w2 * chord / oper.mu

        cl = cl_lerp((aoa, Re))
        cd = cd_lerp((aoa, Re))

        lift = 0.5 * oper.rho * w2**2 * chord * loc.dr * cl
        drag = 0.5 * oper.rho * w2**2 * chord * loc.dr * cd

        thrust = -prop.geometry.n_blades * (lift * np.sin(beta_2) + drag * np.cos(beta_2))
        torque =  prop.geometry.n_blades * (lift * np.cos(beta_2) - drag * np.sin(beta_2)) * loc.r

        # New approximations for the absolute velocity components
        mass_flow = 2 * np.pi * loc.r * loc.dr * oper.rho * va2
        va3_new = oper.v_inf + thrust / mass_flow
        vu2p_new = torque / (mass_flow * loc.r)

        # Stop the iterations if the speeds have sufficiently converged
        va3_ok = np.allclose(va3, va3_new, atol=atol, rtol=rtol)
        vu2p_ok = np.allclose(vu2p, vu2p_new, atol=atol, rtol=rtol)
        if va3_ok and vu2p_ok:
            flag_converged = True
            break

        # Otherwise update the speeds for the next iteration
        vu2p = (1-relax) * vu2p + relax * vu2p_new
        va3 = (1-relax) * va3 + relax * va3_new

    if not flag_converged:
        print(
            f"bemt -- No local convergence after {n_iter} iterations."
            f"\n\t{prop.pretty_name=} -- {loc.r=:.4f} -- {oper.Om=:.0f}"
        )

    return LocalBemSolution(
        prop=prop, oper=oper, loc=loc,
        va2=va2, vu2=vu2, wu2=wu2, thrust=thrust, torque=torque,
        flag_converged=flag_converged
    )


def bem(prop: Propeller, oper: OperatingConditions, sdiv=20) -> BemSolution:
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

    local_bem_sol_list = np.zeros(r_dist.shape, dtype=LocalBemSolution)

    for i, (r, dr) in enumerate(zip(r_dist, dr_dist)):
        local_bem_sol_list[i] = local_bem(prop, oper, Location(r, dr))

    return BemSolution(lsols=local_bem_sol_list)
