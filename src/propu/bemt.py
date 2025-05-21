"""Blade element momentum techniques."""

from typing import NamedTuple

import numpy as np
from numpy.typing import NDArray
from scipy import interpolate

from propu.constant import uconv


class Airfoil(NamedTuple):
    """Tabulated Airfoil data.

    Contain the tabulated coordinates of the airfoil profile,
    as well as the tabulated values of cl and cd,
    for the corresponding sampled angles of attack and Reynolds numbers.
    """

    profile: NDArray[np.float64]  # Sample of airfoil coordinates [m]
    aoa: NDArray[np.float64]  # Sample of angles of attack [rad]
    Re: NDArray[np.float64]  # Sample of Reynolds numbers [-]
    cl: NDArray[np.float64]  # Corresponding lift coefficients [-]
    cd: NDArray[np.float64]  # Corresponding drag coefficients [-]


class PropellerGeometry(NamedTuple):
    span: float  # Span of one blade [m]
    n_blades: int  # Number of blades [-]
    stations: NDArray[np.float64]  # Station locations [m]
    chords: NDArray[np.float64]  # Chords evaluated at stations [m]
    pitches: NDArray[np.float64]  # Pitches evaluated at stations [m]


class Propeller(NamedTuple):
    """Specifications of a generic propeller.

    Attributes
    ----------
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
    """
    keyword: str
    pretty_name: str
    airfoil: Airfoil
    geometry: PropellerGeometry


class OperatingConditions(NamedTuple):
    Om: float  # Rotation speed [rad/s]
    v_inf: float  # Wind speed [m/s]
    rho: float  # Density of the air [kg/m³]
    mu: float  # Dynamic viscosity of the air [Pa*s]


class Location(NamedTuple):
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
    beta: float  # Local relative flow angle [rad]
    thrust: float  # Thrust generated in the local stream tube [N]
    torque: float  # Torque generated in the local stream tube [N*m]
    # Status
    converged: bool  # Raised if the solution has converged.


class BemSolution:
    def __init__(self, lsols: NDArray[LocalBemSolution]):
        # Context
        self.prop = lsols[0].prop  # Which propeller is used ?
        self.oper = lsols[0].oper  # Under what operating conditions ?
        # Extract data from local bem solutions in separate lists
        self.r_dist = np.array([sol.loc.r for sol in lsols])
        self.dr_dist = np.array([sol.loc.dr for sol in lsols])
        self.converged_dist = np.array([sol.converged for sol in lsols])
        self.va2_dist = np.array([sol.va2 for sol in lsols])
        self.vu2_dist = np.array([sol.vu2 for sol in lsols])
        self.wu2_dist = np.array([sol.wu2 for sol in lsols])
        self.beta_dist = np.array([sol.beta for sol in lsols])
        self.thrust_dist = np.array([sol.thrust for sol in lsols])
        self.torque_dist = np.array([sol.torque for sol in lsols])
        self.power_dist = self.torque_dist * self.oper.Om
        # Build global propeller performances quantities
        self.thrust = np.sum(self.thrust_dist)  # Total thrust provided by the propeller [N]
        self.torque = np.sum(self.torque_dist)  # Total torque provided by the propeller [N*m]
        self.power = np.sum(self.power_dist)  # Total power neeed to drive the propeller [W]
        self.eta = self.thrust / self.power * self.oper.v_inf  # Propulsive efficiency [-]


def local_bem(
    prop: Propeller, oper: OperatingConditions, loc: Location,
    *, relax=0.3, atol=1e-4, rtol=1e-3, max_iter=100,
) -> LocalBemSolution:
    """Momentum balance in one stream tube.

    Parameters
    ----------
    prop: Propeller
        The propeller on which the momentum balance is performed.
    oper: OperatingConditions
        The conditions at which the propeller is operating.
    loc: Location
        The location of the blade element (and associated stream tube) on which the local momentum
        balance is performed.
    relax : float, optional
        Relaxation factor that stabilizes the convergence, by reducing the amplitude of the step
        corrections.
    atol, rtol: float, optional
        The maximal absolute and relative difference of the downstream speeds between each steps
        that should be obtained in order to claim a convergence of the solution.
    max_iter: int, optional
        The maximum number of iterations that are tolerated. A warning is printed if no convergence
        is obtained after this number of iterations.
    """
    # Linear interpolation on the blade span to retrieve pitch and chord at `loc`
    chord_lerp = interpolate.make_interp_spline(prop.geometry.stations, prop.geometry.chords, k=1)
    pitch_lerp = interpolate.make_interp_spline(prop.geometry.stations, prop.geometry.pitches, k=1)
    chord = chord_lerp(loc.r)
    pitch = pitch_lerp(loc.r)

    # Linear interpolator methods for lift and drag coefficients
    cl_lerp = interpolate.RegularGridInterpolator(
        (prop.airfoil.aoa, prop.airfoil.Re),
        prop.airfoil.cl,
        bounds_error=False,
        fill_value=None,
    )
    cd_lerp = interpolate.RegularGridInterpolator(
        (prop.airfoil.aoa, prop.airfoil.Re),
        prop.airfoil.cd,
        bounds_error=False,
        fill_value=None,
    )

    # Initial guesstimates and iteration parameters
    #
    # NOTE: Avoid to set null speed for va3
    # If this is the case, then the computed mass flow in null,
    # hence va3_new and vu2p_new are NaN.
    va3 = max(oper.v_inf, 1.0)
    vu2p = 0.0  # [m/s]
    converged = False  # Raised if the speeds have been stabilized

    # Determine forces and velocities iteratively
    for n_iter in range(1, max_iter + 1):
        # Speed triangle at the blade
        va2 = (oper.v_inf + va3) / 2
        vu2 = vu2p / 2
        wu2 = vu2 - oper.Om * loc.r
        w2 = np.sqrt(va2**2 + wu2**2)
        beta = np.atan2(wu2, va2)

        # Lift and drag of one blade
        aoa = pitch - (np.pi / 2 + beta)
        Re = oper.rho * w2 * chord / oper.mu
        cl = cl_lerp((aoa, Re))
        cd = cd_lerp((aoa, Re))
        lift = 0.5 * oper.rho * w2**2 * chord * loc.dr * cl
        drag = 0.5 * oper.rho * w2**2 * chord * loc.dr * cd

        # Momentum balance on the streamtube,
        # leading to new estimates of the absolute velocity components
        thrust = -prop.geometry.n_blades * (lift * np.sin(beta) + drag * np.cos(beta))
        torque = prop.geometry.n_blades * (lift * np.cos(beta) - drag * np.sin(beta)) * loc.r
        mass_flow = 2 * np.pi * loc.r * loc.dr * oper.rho * va2
        va3_new = oper.v_inf + thrust / mass_flow
        vu2p_new = torque / (mass_flow * loc.r)

        # Stop the iterations if the speeds have sufficiently converged
        va3_ok = np.allclose(va3, va3_new, atol=atol, rtol=rtol)
        vu2p_ok = np.allclose(vu2p, vu2p_new, atol=atol, rtol=rtol)
        if va3_ok and vu2p_ok:
            converged = True
            break

        # Otherwise update the speeds for the next iteration
        vu2p = (1 - relax) * vu2p + relax * vu2p_new
        va3 = (1 - relax) * va3 + relax * va3_new

    if not converged:
        print(
            f"bemt -- No local convergence after {n_iter} iterations."
            f"\n\tpropeller: {prop.pretty_name}"
            f" -- radius: {loc.r:.4f} m"
            f" -- rotation: {oper.Om * uconv('rad/s', 'rpm'):.0f} rpm"
        )

    return LocalBemSolution(
        prop=prop, oper=oper, loc=loc,
        va2=va2, vu2=vu2, wu2=wu2, beta=beta,
        thrust=thrust, torque=torque,
        converged=converged
    )


def bem(prop: Propeller, oper: OperatingConditions, sdiv=20, **local_bem_kwargs) -> BemSolution:
    """Calculate quantities derived from momentum balance of a propeller.

    This function applies the blade element method to the propeller `prop`, for the specified
    operating conditions `oper`, and the number of subdivisions `sdiv` in the blades.
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
        local_bem_sol_list[i] = local_bem(prop, oper, Location(r, dr), **local_bem_kwargs)

    return BemSolution(lsols=local_bem_sol_list)
