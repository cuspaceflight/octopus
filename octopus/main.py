"""Implementation of injector fluid dynamics classes.

"""
from dataclasses import dataclass
from typing import Iterable, Sequence

import CoolProp.CoolProp as CP
import numpy as np
from CoolProp import AbstractState
from CoolProp.CoolProp import PropsSI

from .utils import fd


class Fluid:
    """Represents a fluid."""

    def __init__(self, name: str):
        """Initiate a Fluid instance.

        :param ID: id of fluid
        """

        self.name = name
        self.state = AbstractState('HEOS', name)
        self.Tmax = self.state.Tmax()
        self.Tmin = self.state.Tmin()
        self.pmax = self.state.pmax()
        self.pmin = 0

    def copy(self):
        fluid = Fluid(self.name)
        fluid.set_state(P=self.state.p(), T=self.state.T())
        return fluid

    def set_state(self, D=None, P=None, T=None, Q=None, H=None, S=None, U=None):
        if sum([bool(D), bool(P), bool(T), bool(Q), bool(H), bool(S), bool(U)]) != 2:
            raise ValueError('Must have exactly 2 arguments')
        args = []
        if D:
            if P:
                args = CP.DmassP_INPUTS, D, P
            elif T:
                args = CP.DmassT_INPUTS, D, T
            elif Q:
                args = CP.DmassQ_INPUTS, D, Q
            elif H:
                args = CP.DmassHmass_INPUTS, D, H
            elif S:
                args = CP.DmassSmass_INPUTS, D, S
            elif U:
                args = CP.DmassUmass_INPUTS, D, U
        elif P:
            if T:
                args = CP.PT_INPUTS, P, T
            elif Q:
                args = CP.PQ_INPUTS, P, Q
            elif H:
                args = CP.HmassP_INPUTS, H, P
            elif S:
                args = CP.PSmass_INPUTS, P, S
            elif U:
                args = CP.PUmass_INPUTS, P, U
        elif T:
            if Q:
                args = CP.QT_INPUTS, Q, T
            elif H:
                args = CP.HmassT_INPUTS, H, T
            elif S:
                args = CP.SmassT_INPUTS, S, T
            elif U:
                args = CP.TUmass_INPUTS, T, U
        elif Q:
            if H:
                args = CP.HmassQ_INPUTS, H, Q
            elif S:
                args = CP.QSmass_INPUTS, Q, S
            elif U:
                raise ValueError('Invalid combination: Q and U')
        elif H:
            if S:
                args = CP.HmassSmass_INPUTS, H, S
            elif U:
                raise ValueError('Invalid combination: H and U')
        elif S and U:
            args = CP.SmassUmass_INPUTS, S, U
        else:
            raise ValueError('Invalid combination')

        try:
            self.state.update(*args)
        except ValueError:
            pass

    def Vf(self, S: float):
        x = self.state.Q()
        if x == -1:
            if self.state.rhomass() > self.rhol([self.state.T()])[0]:
                return 0
            else:
                return 1
        else:
            return 1 / (1 / (1 + S * ((1 - x) / x) * self.rhog([self.state.T()])[0] / self.rhol([self.state.T()])[0]))

    def chi(self):
        x = self.state.Q()
        if x == -1:
            if self.state.rhomass() > self.rhol([self.state.T()])[0]:
                return 0
            else:
                return 1
        else:
            return x

    def viscosityl(self, T):
        b1 = 1.6089
        b2 = 2.0439
        b3 = 5.24
        b4 = 0.0293423
        theta = (self.state.T_critical() - b3) / (T - b3)
        return b4 * np.exp(b1 * (theta - 1) ** (1 / 3) + b2 * (theta - 1) ** (4 / 3))

    def viscosityg(self, T):
        b1 = 3.3281
        b2 = -1.18237
        b3 = -0.055155
        Tr = T / self.state.T_critical()
        return 0.001 * np.exp(b1 + b2 * (1 / Tr - 1) ** (1 / 3) + b3 * (1 / Tr - 1) ** (4 / 3)) if (
                self.Tmin < T < self.Tmax) else None

    def viscosity(self):
        a = self.Vf(S=1)
        T = self.state.T()
        if a == 0:
            return self.viscosityl(T)
        if a == 1:
            return self.viscosityg(T)
        mug = self.viscosityg(T)
        mul = self.viscosityl(T)
        x = self.state.Q()
        rho = self.state.rhomass()
        rhol = self.rhol([T])[0]
        rhog = self.rhog([T])[0]
        return mug * x * rho / rhog + mul * (1 - x) * rho / rhol

    def dvg_dP_saturation(self):
        dP = 0.002
        rho1 = PropsSI('D', 'P', self.state.p() + dP / 2, 'Q', 1, self.name)
        rho0 = PropsSI('D', 'P', self.state.p() - dP / 2, 'Q', 1, self.name)
        dvg = 1 / rho1 - 1 / rho0
        return dvg / dP

    def rhol(self, T: Iterable):
        return [PropsSI('D', 'T', t, 'Q', 0, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def rhog(self, T: Iterable):
        return [PropsSI('D', 'T', t, 'Q', 1, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def psat(self, T: Iterable):
        return [PropsSI('P', 'Q', 0.5, 'T', t, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def tsat(self, P: Iterable):
        return [PropsSI('T', 'P', P, 'Q', 0.5, self.name) if (self.pmin < p < self.pmax) else None for p in P]

    def hl(self, T: Iterable):
        return [PropsSI('H', 'T', t, 'Q', 0, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def hg(self, T: Iterable):
        return [PropsSI('H', 'T', t, 'Q', 1, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def sl(self, T: Iterable):
        return [PropsSI('S', 'T', t, 'Q', 0, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def sg(self, T: Iterable):
        return [PropsSI('S', 'T', t, 'Q', 1, self.name) if (self.Tmin < t < self.Tmax) else None for t in T]

    def __repr__(self):
        return f'Fluid({self.name}: p={self.state.p() / 100000:.1f}bar, t={self.state.T():.1f}K)'


@dataclass
class Cell:
    p: float
    x: float


class Mesh:
    def __init__(self, fluid: Fluid, velocity: float, L: float, n: int):
        self.cells = [Cell(fluid=fluid.copy(), velocity=velocity, dx=L / n) for _ in range(n)]


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

    def __init__(self, manifold: Manifold, L: float, D: float = None, A: float = None, orifice_type: int = 0,
                 Cd: float = 0.7):
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

        if not A and not D:
            raise TypeError
        elif not D:
            self.D = np.sqrt(4 * A / np.pi)
            self.A = A
        else:
            self.A = 0.25 * np.pi * D ** 2
            self.D = D

        self.Cd = Cd

        # defines default Orifice.m_dot(p1) function
        self.m_dot = self.m_dot_dyer

    def set_A(self, A):
        self.A = A

    # @lru_cache(maxsize=1)
    def m_dot_SPI(self, p1: float):
        """Return single-phase-incompressible mass flow rate.

        :param p1: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        if self.orifice_type == self.STRAIGHT:
            # find initial conditions
            # T0,p0 known
            p0 = self.manifold.p
            T0 = self.manifold.T
            fluid = self.manifold.fluid

            # check pressure drop
            if p0 < p1:
                return None
            fluid.set_state(T=T0, P=p0)
            rho0 = fluid.state.rhomass()

            # compute mass flow
            return self.A * self.Cd * np.sqrt(2 * rho0 * (p0 - p1))
        else:
            raise NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    # @lru_cache(maxsize=1)
    def m_dot_HEM(self, p1: float):
        """Return homogeneous-equilibrium-model mass flow rate.

        :param p1: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        if self.orifice_type == self.STRAIGHT:

            # find initial conditions
            # p0,T0 known
            p0 = self.manifold.p
            T0 = self.manifold.T
            fluid = self.manifold.fluid

            # check for pressure drop
            if p0 < p1:
                return None

            # retrieve enthalpy and entropy
            fluid.set_state(T=T0, P=p0)
            h0 = fluid.state.hmass()
            s0 = fluid.state.smass()

            # set final conditions
            # p1, s1 known
            p1 = p1
            s1 = s0

            # compute isentropic flash expansion
            fluid.set_state(P=p1, S=s1)
            rho1 = fluid.state.rhomass()
            h1 = fluid.state.hmass()

            # check decrease in enthalpy
            if h1 > h0:
                return None

            # compute mass flow
            return self.A * self.Cd * rho1 * np.sqrt(2 * (h0 - h1))

        else:
            return NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    # @lru_cache(maxsize=1)
    def m_dot_dyer(self, P_cc: float):
        """Return Dyer model mass flow rate.

        :param p0: optio
        :param P_cc: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """

        p0 = self.manifold.p
        T0 = self.manifold.T
        p1 = P_cc

        pv = self.manifold.fluid.psat([T0])[0]

        kappa = np.sqrt((p0 - p1) / (pv - p1))
        W = 1 / (1 + kappa)
        if not (self.m_dot_SPI(P_cc) and self.m_dot_HEM(P_cc)):
            return None
        return (1 - W) * self.m_dot_SPI(P_cc) + W * self.m_dot_HEM(P_cc)

    def dP_HEM_frictional(self, mdot: float):
        # find initial conditions
        # p0,T0 known
        p0 = self.manifold.p
        T0 = self.manifold.T
        fluid = self.manifold.fluid

        # retrieve void fraction
        fluid.set_state(T=T0, P=p0)
        s0 = fluid.state.smass()
        num_cells = 100
        dz = self.L / num_cells
        G = mdot / self.A

        rho0 = fluid.state.rhomass()
        V0 = G/rho0
        dP0 = 0.5*rho0*V0**2

        cells = [Cell(p=p0-dP0, x=0) for _ in range(num_cells)]
        cells = [Cell(p=p0-dP0, x=0), *cells]

        try:

            for i in range(num_cells):
                vf = 1 / fluid.rhol([fluid.state.T()])[0]
                vfg = 1 / fluid.rhog([fluid.state.T()])[0] - vf

                fluid.set_state(P=cells[i].p, S=s0)
                cells[i + 1].x = fluid.chi()

                dx_dz = (cells[i + 1].x - cells[i].x) / dz
                dvg_dP = fluid.dvg_dP_saturation()

                V = G / fluid.state.rhomass()
                Re = fluid.state.rhomass() * V * self.D / fluid.viscosity()

                dP_dz = -((2 * (fd(Re) / 4) * G ** 2 * vf / self.D) * (1 + cells[i + 1].x * (vfg / vf)) + G ** 2 * vf * (
                        vfg / vf) * dx_dz) / (1 + G ** 2 * cells[i + 1].x * dvg_dP)

                cells[i + 1].p = cells[i].p + dP_dz * dz
                print(f'{dP_dz * dz:7.0f}\t{Re:4.0f}\t{fluid.viscosity():2.3f}\t{vf:.5f}\t{vfg:.5f}')

        except Exception:
            pass

        print(dP0,cells[-1].p)
        return [cell.p for cell in cells]


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
