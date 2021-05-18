"""Implementation of injector fluid dynamics classes.

"""
from functools import lru_cache
from typing import Sequence, Iterable

import numpy as np
from CoolProp import AbstractState
from CoolProp.CoolProp import PropsSI


class Fluid:
    """Represents a fluid."""

    def __init__(self, name: str):
        """Initiate a Fluid instance.

        :param ID: id of fluid
        """

        self.name = name
        self.state = AbstractState('HEOS', name)

    def rhol(self, T: Iterable):
        return [PropsSI('D', 'T', t, 'Q', 0, self.name) if (182.23 < t < 309.52) else None for t in T]

    def rhog(self, T: Iterable):
        return [PropsSI('D', 'T', t, 'Q', 1, self.name) if (182.23 < t < 309.52) else None for t in T]

    def psat(self,T:Iterable):
        return [PropsSI('P', 'Q', 0.5, 'T', t, self.name) if (182.23 < t < 309.52) else None for t in T]


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

    def __init__(self, fluid: str, parent: PropertySource):
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
        self.A = 0.2 * np.pi * self.D ** 2
        self.Cd = Cd

        # defines default Orifice.m_dot(p1) function
        self.m_dot = self.m_dot_dyer

    @lru_cache(maxsize=1)
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

            rho0 = PropsSI('D', 'T', T0, 'P', p0, fluid)

            # compute mass flow
            return self.A * np.sqrt(2 * rho0 * (p0 - p1))
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
            # p0,T0 known
            p0 = self.manifold.p
            T0 = self.manifold.T
            fluid = self.manifold.fluid

            # check for pressure drop
            if p0 < p1:
                return None

            # retrieve enthalpy and entropy
            h0 = PropsSI('H', 'T', T0, 'P', p0, fluid)
            s0 = PropsSI('S', 'T', T0, 'P', p0, fluid)

            # set final conditions
            # p1, s1 known
            p1 = p1
            s1 = s0

            # compute isentropic flash expansion
            rho1 = PropsSI('D', 'P', p1, 'S', s1, fluid)
            h1 = PropsSI('H', 'P', p1, 'S', s1, fluid)

            # check decrease in enthalpy
            if h1 > h0:
                return None

            # compute mass flow
            return self.A * rho1 * np.sqrt(h0 - h1)

        else:
            return NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    @lru_cache(maxsize=1)
    def m_dot_dyer(self, P_cc: float):
        """Return Dyer model mass flow rate.

        :param P_cc: combustion chamber pressure (Pa)
        :return: mass flow rate (kg/s)

        """
        p0 = self.manifold.p
        T0 = self.manifold.T
        p1 = P_cc
        fluid = self.manifold.fluid

        pv = PropsSI('P', 'Q', 0.5, 'T', T0, fluid)

        kappa = np.sqrt((p0 - p1) / (pv - p1))
        W = 1 / (1 + kappa)
        if not (self.m_dot_SPI(P_cc) and self.m_dot_HEM(P_cc)):
            return None
        return self.Cd * ((1 - W) * self.m_dot_SPI(P_cc) + W * self.m_dot_HEM(P_cc))


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
