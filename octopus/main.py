"""Implementation of injector fluid dynamics classes.

References:
    - [1]   Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
    - [2]   An investigation of injectors for use with high vapor pressure propellants with applications to hybrid rockets,
            Waxman B.S., https://stacks.stanford.edu/file/druid:ng346xh6244/BenjaminWaxmanFinal-augmented.pdf
    - [3]   Short Fundamental Equations of State for 20 Industrial Fluids,
            Lemmon; Span, https://pubs.acs.org/doi/pdf/10.1021/je050186n
    - [4]   Engineering Model to Calculate Mass Flow Rate of a Two-Phase Saturated Fluid Through An Injector Orifice,
            Solomon B.J., https://digitalcommons.usu.edu/cgi/viewcontent.cgi?article=1110&context=gradreports

"""
from functools import lru_cache
from json import load
from os.path import dirname
from typing import List, Dict, Sequence

import thermo.chemical
from numpy import pi, sqrt, array, log, exp, nan_to_num, inf
from scipy.optimize import least_squares

from .utils import derivative


class Fluid:
    """Represents a fluid with a choice of Helmholz or thermo EOS methods"""

    def __init__(self, ID: str, T: float = 298.15, P: float = 101325, method='thermo'):
        """Initiate a Fluid instance.

        :param ID: id of fluid with a name recognisable by thermo
        :param T: Temperature (K)
        :param P: Pressure (Pa)
        """
        chem = thermo.chemical.Chemical(ID, T, P)
        self.method = method.lower()

        if self.method == 'helmholz':
            with open(f'{dirname(__file__)}/data/{chem.CAS}.json', 'r') as f:
                data = load(f)

            self.a = data['a']
            self.c = data['c']
            self.v_ = array(data['v'])
            self.u_ = array(data['u'])
            self.n = data['n']
            self.Ag = data['Ag']
            self.Al = data['Al']
            self.polar = bool(data['polar'])
            self.Tc = data['Tc']
            self.rhoc = data['rhoc']
            self.R_specific = data['R_specific']
            self.VaporPressure = chem.VaporPressure

        elif self.method == 'thermo':
            pass
        else:
            raise NotImplementedError(f'method cannot be: {self.method}')

    @lru_cache(maxsize=1)
    def alpha_0(self, delta: float, tau: float):
        """Calculate the reduced ideal gas Helmholz energy.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: ideal gas component of Helmholz energy
        :rtype: float
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        alpha_0 = (self.a[0]
                   + self.a[1] * tau
                   + log(delta)
                   + (self.c[0] - 1) * log(tau)
                   - (self.c[1] * self.Tc ** self.c[2]) / (self.c[2] * (self.c[2] + 1)) * tau ** (-self.c[2])
                   + sum(self.v_ * log(1 - exp(-self.u_ * tau / self.Tc))))
        return alpha_0

    @lru_cache(maxsize=1)
    def alpha_r(self, delta: float, tau: float) -> float:
        """Calculate the reduced helmholz residual energy.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: residual component of Helmholz energy
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        if self.polar:

            alpha_r = (self.n[0] * delta * tau ** 0.25
                       + self.n[1] * delta * tau ** 1.25
                       + self.n[2] * delta * tau ** 1.5
                       + self.n[3] * delta ** 3 * tau ** 0.25
                       + self.n[4] * delta ** 7 * tau ** 0.875
                       + self.n[5] * delta * tau ** 2.375 * exp(-delta)
                       + self.n[6] * delta ** 2 * tau ** 2 * exp(-delta)
                       + self.n[7] * delta ** 5 * tau ** 2.125 * exp(-delta)
                       + self.n[8] * delta * tau * 3.5 * exp(-(delta ** 2))
                       + self.n[9] * delta * tau ** 6.5 * exp(-(delta ** 2))
                       + self.n[10] * delta ** 4 * tau ** 4.75 * exp(-(delta ** 2))
                       + self.n[11] * delta ** 2 * tau ** 12.5 * exp(-(delta ** 3)))
            return alpha_r
        else:
            raise NotImplementedError

    def rho_g(self, T: float) -> float:
        """Evaluate vapour density at T.

        :param T: Temperature (K)
        :return: vapour density
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        tau = 1 - T / self.Tc
        return nan_to_num(self.rhoc + self.Ag[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(self.Ag)))

    def rho_l(self, T: float) -> float:
        """Evaluate liquid density at T.

        :param T: Temperature (K)
        :return: liquid density
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        tau = 1 - T / self.Tc
        return nan_to_num(self.rhoc + self.Al[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(self.Al)))

    @lru_cache(maxsize=1)
    def a0_t(self, delta: float, tau: float) -> float:
        """Evaluate derivative of alpha_0 with respect to tau.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d(alpha_0)/d(tau)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.alpha_0, 1, delta, tau)

    @lru_cache(maxsize=1)
    def a0_tt(self, delta: float, tau: float) -> float:
        """Evaluate second derivative of alpha_0 with respect to tau.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d2(alpha_0)/d(tau)d(tau)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.a0_t, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_d(self, delta: float, tau: float) -> float:
        """Evaluate derivative of alpha_r with respect to delta.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d(alpha_r)/d(delta)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.alpha_r, 0, delta, tau)

    @lru_cache(maxsize=1)
    def ar_t(self, delta: float, tau: float) -> float:
        """Evaluate derivative of alpha_r with respect to tau.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d(alpha_r)/d(tau)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.alpha_r, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_dd(self, delta: float, tau: float) -> float:
        """Evaluate second derivative of alpha_r with respect to delta.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d2(alpha_r)/d(delta)d(delta
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.ar_d, 0, delta, tau)

    @lru_cache(maxsize=1)
    def ar_dt(self, delta: float, tau: float) -> float:
        """Evaluate derivative of alpha_r with respect to delta and tau.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d2(alpha_r)/d(delta)d(tau)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.ar_d, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_tt(self, delta: float, tau: float) -> float:
        """Evaluate second derivative of alpha_r with respect to tau.

        :param delta: rho/rho_c to be evaluated at
        :param tau: T_c/T to be evaluated at
        :return: d2(alpha_r)/d(tau)d(tau)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        return derivative(self.ar_t, 1, delta, tau)

    @lru_cache(maxsize=1)
    def dp_drho(self, rho: float, T: float) -> float:
        """Evaluate derviative of pressure with respect to density.

        :param rho: density to be evaluated at
        :param T: temperature to be evaluated at
        :return: d(p)/d(rho)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * (1 + 2 * delta * self.ar_d(delta, tau) + delta * delta * self.ar_dd(delta, tau))

    @lru_cache(maxsize=1)
    def dg_drho(self, rho: float, T: float) -> float:
        """Evaluate derviative of specific gibbs free energy with respect to density.

        :param rho: density to be evaluated at
        :param T: temperature to be evaluated at
        :return: d(g)/d(rho)
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * (2 * self.ar_d(delta, tau) + delta * self.ar_dd(delta, tau))

    @lru_cache(maxsize=1)
    def p(self, rho: float, T: float) -> float:
        """Evaluate pressure.

        :param rho: density
        :param T: temperature
        :return: pressure
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return rho * self.R_specific * T * (1 + delta * self.ar_d(delta, tau))

    @lru_cache(maxsize=1)
    def u(self, rho: float, T: float) -> float:
        """Evaluate specific internal energy.

        :param rho: density
        :param T: temperature
        :return: specific internal energy
        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau))

    @lru_cache(maxsize=1)
    def h(self, rho: float, T: float) -> float:
        """Evaluate specific enthalpy.

        :param rho: density
        :param T: temperature
        :return: specific enthalpy

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return (self.R_specific * T
                * (1 + tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) + delta * self.ar_d(delta, tau)))

    @lru_cache(maxsize=1)
    def s(self, rho: float, T: float) -> float:
        """Evaluate specific entropy.

        :param rho: density
        :param T: temperature
        :return: specific entropy

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return (self.R_specific *
                (tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) - self.alpha_0(delta, tau) - self.alpha_r(delta,
                                                                                                                 tau)))

    @lru_cache(maxsize=1)
    def g(self, rho: float, T: float) -> float:
        """Evaluate specific gibbs free energy

        :param rho: density
        :param T: temperature
        :return: specific gibbs free energy

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return (self.R_specific * T
                * (1 + self.alpha_0(delta, tau) + self.alpha_r(delta, tau) + delta * self.ar_d(delta, tau)))

    @lru_cache(maxsize=1)
    def z(self, rho: float, T: float) -> float:
        """Return compressibility

        :param rho: density
        :param T: temperature
        :return: compressibility, z

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')
        delta = rho / self.rhoc
        tau = self.Tc / T
        return 1 + delta * self.ar_d(delta, tau)

    @lru_cache(maxsize=1)
    def cp(self, rho: float, T: float) -> float:
        """Return specific heat capacity at constant pressure.

        :param rho: density
        :param T: temperature
        :return: heat capacity at constant pressure

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')

        return derivative(self.h, 1, rho, T)  # cp = dh/dT

    @lru_cache(maxsize=1)
    def cv(self, rho: float, T: float) -> float:
        """Return specific heat capacity at constant volume.

        :param rho: density
        :param T: temperature
        :return: heat capacity at constant volume

        """
        if self.method != 'helmholz':
            raise NotImplementedError(f'Helmholz energy methods not available for EOS: {self.method}')

        return derivative(self.u, 1, rho, T)  # cv = du/dT

    def fun_ps(self, x: Sequence[float], u: Sequence[str], y: Sequence[float]) -> List[float]:
        """Return a vectorised cost function to be used in a property solver.

        See Usage for more info.

        :param x: object containing arguments to be passed to get_properties
        :param u: object containing property names of properties to be returned and compared
        :param y: object containing property values to compare returned values to
        :return: object containing difference between calculated property values and those in y

        """
        if len(x) != len(u) != len(y) != 2:
            raise ValueError('can only solve for 2 variables, from 2 knowns')
        if 'rho' in u and 'T' in u:
            return [0, 0]

        if 'rho' in u:
            i = u.index('rho')
            var = u[i - 1]
            val = y[i - 1]
            res = [self.get_properties(x[0], x[1])[var] - val]
            res.insert(i, x[0] - y[i])
            return res

        elif 'T' in u:
            i = u.index('T')
            var = u[i - 1]
            val = y[i - 1]
            res = [self.get_properties(x[0], x[1])[var] - val]
            res.insert(i, x[1] - y[i])
            return res

        else:
            return [self.get_properties(x[0], x[1])[var] - val for var, val in zip(u, y)]

    @lru_cache(maxsize=1)
    def get_properties(self, rho: float, T: float) -> Dict[str, float]:
        """Return properties in a dict

        :param rho: density
        :param T: temperature
        :return: dict where keys are property names referring to float values
        """
        rhog = self.rho_g(T)
        rhol = self.rho_l(T)

        if rho < rhog:
            p = self.p(rho, T)
            chi = 1.0
            h = self.h(rho, T)
            s = self.s(rho, T)

        elif rho > rhol:
            p = self.p(rho, T)
            chi = 0.0
            h = self.h(rho, T)
            s = self.s(rho, T)
        else:
            p = self.VaporPressure(T)
            chi = rhog * (rhol - rho) / (rho * (rhol - rhog))
            h = self.h(rhog, T) * chi + self.h(rhol, T) * (1 - chi)
            s = self.s(rhog, T) * chi + self.s(rhol, T) * (1 - chi)

        res = {'p': p, 'chi': chi, 'h': h, 's': s}

        return res


class PropertySource:
    """An object to provide a constant pressure and temperature as a parent object for a :class:`Manifold`.

    The user may create their own property source class, by extending :class`PropertySource`, and modifying the current
    :meth:`p` and :meth:`T` methods."""

    def __init__(self, p: float = 101325, T: float = 298):
        """Initialise :class`PropertySource` object with constant pressure and temperature to supply to manifold fluid.
        :param p: pressure (Pa)
        :param T: temperature (K)
        """
        self._p = p
        self._T = T

    @property
    def p(self):
        return self._p

    @property
    def T(self):
        return self._T

    def __repr__(self):
        print(f'PropertySource: p = {self._p}, T = {self._T}')


class Manifold:
    """Represent a propellant manifold, at least one :class:`Fluid` input and one
    :class:`Orifice` output. If a user wishes to model losses within the manifold, they may extend this
    class and add the required computing into the :meth:`Manifold.p` and :meth:`Manifold.T`
    functions. """

    def __init__(self, fluid: Fluid, parent: PropertySource):
        """Initialise :class:`Manifold` object with a working fluid and property parent.

        :param fluid: :class:`Fluid` object to use EOS functions from
        :param parent: :class:`PropertySource` object to get p and T from
        """
        self.fluid = fluid
        self.parent = parent
        self.chi = 0

    @property
    def p(self):
        return self.parent.p

    @property
    def T(self):
        return self.parent.T


class Orifice:
    """Represent a propellant orifice on the injector plate, at least one :class:`Manifold` input."""

    STRAIGHT = 0
    CAVITATING = 1

    def __init__(self, manifold: Manifold, L: float, D: float, orifice_type: int = 0, Cd: float = 0.7):
        """Initialise :class:`Orifice` object.

        :param manifold: :class:`Manifold` to get fluid EOS and properties from
        :param L: orifice length (m)
        :param D: orifice diameter (m)
        :param orifice_type: :attr:`STRAIGHT` (default) or :attr:`CAVITATING`
        :param Cd: discharge coefficient = 0.7

        """
        # manifold holds fluid info
        self.manifold = manifold

        self.orifice_type = orifice_type
        self.L = L
        self.D = D
        self.A = 0.2 * pi * self.D ** 2
        self.Cd = Cd

        # defines default Orifice.m_dot(p1) function
        if manifold.fluid.method == 'helmholz':
            self.m_dot = self.m_dot_dyer
        else:
            self.m_dot = self.m_dot_SPI

    @lru_cache(maxsize=1)
    def m_dot_SPI(self, p1: float):
        """Return single-phase-incompressible mass flow rate.

        :param p1: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        if self.orifice_type == self.STRAIGHT:
            # find initial conditions
            # chi0,p0 known
            p0 = self.manifold.p
            chi0 = self.manifold.chi
            fluid = self.manifold.fluid

            # check pressure drop
            if p0 < p1:
                return None

            # use helmholz EOS if available
            if fluid.method == 'helmholz':
                u = ['p', 'chi']
                y = [p0, chi0]

                # compute x resulting in y
                initial = least_squares(fluid.fun_ps, [800, 250], bounds=([0, 0], [inf, fluid.Tc]), args=[u, y])
                rho0, T0 = initial.x

            elif fluid.method == 'thermo':
                # TODO: integrate rho functions
                rho0 = fluid.rhol
            else:
                raise NotImplementedError(f'EOS method "{fluid.method}" not implemented')

            # compute mass flow
            return self.A * sqrt(2 * rho0 * (p0 - p1))
        else:
            raise NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    @lru_cache(maxsize=1)
    def m_dot_HEM(self, p1: float):
        """Return homogeneous-equilibrium-model mass flow rate.

        :param p1: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        if self.orifice_type == self.STRAIGHT:

            # find initial conditions
            # chi0,p0 known
            p0 = self.manifold.p
            chi0 = self.manifold.chi
            fluid = self.manifold.fluid

            # check for pressure drop
            if p0 < p1:
                return None

            # check fluid has a helmholz eos
            if fluid.method != 'helmholz':
                raise NotImplementedError(f'EOS method must be "helmholz", not "{fluid.method}" for HEM modelling')

            # set known values to iterate towards
            u = ['p', 'chi']
            y = [p0, chi0]

            # compute x resulting in y
            initial = least_squares(fluid.fun_ps, [800, 250], bounds=([0, 0], [inf, fluid.Tc]), args=[u, y])
            rho0, T0 = initial.x

            # retrieve enthalpy and entropy
            h0 = fluid.get_properties(rho0, T0)['h']
            s0 = fluid.get_properties(rho0, T0)['s']

            # set final conditions
            # p1, s1 known
            p1 = p1
            s1 = s0

            # set known values to iterate towards
            u = ['p', 's']
            y = [p1, s1]

            # compute x resulting in y
            final = least_squares(fluid.fun_ps, [500, 245], bounds=([0, 0], [inf, fluid.Tc]), args=[u, y])
            rho1, T1 = final.x

            # retrieve enthalpy
            h1 = fluid.get_properties(rho1, T1)['h']

            # check decrease in enthalpy
            if h1 > h0:
                return None

            # compute mass flow
            return self.A * rho1 * sqrt(h0 - h1)

        else:
            return NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    @lru_cache(maxsize=1)
    def m_dot_dyer(self, P_cc: float):
        """Return Dyer model mass flow rate.

        :param P_cc: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        kappa = 1  # fluid enters orifice at vapour pressure
        W = 1 / (1 + kappa)
        if not (self.m_dot_SPI(P_cc) and self.m_dot_HEM(P_cc)):
            return None
        return self.Cd * ((1 - W) * self.m_dot_SPI(P_cc) + W * self.m_dot_HEM(P_cc))

    @lru_cache(maxsize=1)
    def m_dot_waxman(self, p1: float, choke_margin: float = 0.8,
                     Cd_inj: float = 0.65, Cd_diff: float = 0.99, suppress=False):
        """Return estimate of cavitating injector mass flow rate and throat diameter.

        See Ref [2], section 6.

        D is used as the exit diameter.
        Both discharge coefficients default to 1 if not specified.
        Recommend 0.65 for Cd_total, 0.99 for Cd_diff.

        To increase the probability of the onset of choking, the safety
        factor should be increased.

        To be used as a tool for informing cold flow tests, NOT accurate simulation.

        :param p1: combustion chamber pressure (Pa)
        :param choke_margin: rough margin to reduce the vapour pressure why
        :param Cd_inj: injector discharge coefficient
        :param Cd_diff: diffuser discharge coefficient
        :return: mass flow rate (kg/s)

        """
        fluid = self.manifold.fluid
        chi0 = self.manifold.chi

        if self.orifice_type != 1:
            raise NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented for Waxman orifice analysis')

        if fluid.method != 'helmholz':
            raise NotImplementedError(f'EOS method "{fluid.method}" not implemented')

        P1 = self.manifold.p

        P2 = p1
        A2 = self.A

        if p1 > P1:
            raise ValueError("Downstream pressure exceeeds upstream pressure")

        u = ['p', 'chi']
        y = [P1, chi0]

        initial = least_squares(fluid.fun_ps, [800, 250], bounds=([0, 0], [inf, fluid.Tc]), args=[u, y])
        rho1, T1 = initial.x

        Pv = fluid.Psat

        if Pv * choke_margin > p1 and suppress is False:
            print("Cannot choke flow: Chill propellant or increase downstream pressure.")

        # allow for fluid cooling in injector
        Pt_target = Pv * choke_margin

        # Calculate the throat area from Ref [2], equation 6.33
        At = A2 / sqrt(1 + (Cd_diff / Cd_inj) * (P2 - Pt_target) / (P1 - P2))
        Dt = sqrt(4 * At / pi)

        # Calculate actual throat pressure from Ref [2], equation 6.21
        Pt_actual = P2 - (Cd_inj / Cd_diff) ** 2 * (A2 / At) ** 2 * (P1 - P2) * (1 - (At / A2) ** 2)

        # mdot_diff is the choked mass flow through the diffuser section.
        mdot_diff = Cd_diff * At * sqrt(2 * rho1 * (p1 - Pt_actual) / (1 - (At / A2) ** 2))
        mdot_inj = Cd_inj * A2 * sqrt(2 * rho1 * (P1 - p1))

        # Perform sense checks (by continuity, mdot_diff = mdot_inj)
        # Unless intentionally ignored
        if suppress is False:
            if round(mdot_diff, 6) != round(mdot_inj, 6):
                raise Warning("Continuity check failed on Waxman injector.")

            if round(Pt_actual, 6) != round(Pt_target, 6):
                print("Target throat pressure was not acheived in Waxman injector.")

        return {"m_dot": mdot_inj, "Dt": Dt}


class Element:
    """Represent an injector element, at least one :class:`Orifice` input"""

    def __init__(self, o_orifices: Sequence[Orifice], f_orifices: Sequence[Orifice]):
        self.o_orifices = o_orifices
        self.f_orifices = f_orifices

    def of_ratio(self, p: float):
        m_dot_o = 0
        for o in self.o_orifices:
            m = o.m_dot(p)
            if m:
                m_dot_o += m
            else:
                return None
        m_dot_f = 0
        for f in self.f_orifices:
            m = f.m_dot(p)
            if m:
                m_dot_f += m
            else:
                return None
        return m_dot_o / m_dot_f
