"""Implementation of injector fluid dynamics classes.

"""
from functools import lru_cache
from typing import Iterable, Sequence

import CoolProp.CoolProp as CP
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
        self.Tmax = self.state.Tmax()
        self.Tmin = self.state.Tmin()
        self.pmax = self.state.pmax()
        self.pmin = 0

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
            raise AttributeError("need either A or D as an input")
        elif A:
            self.A = A
        else:
            self.A = 0.25 * np.pi * D ** 2

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
            return self.A * self.Cd * rho1 * np.sqrt(2*(h0 - h1))

        else:
            return NotImplementedError(f'Orifice type "{self.orifice_type}" not implemented')

    # @lru_cache(maxsize=1)
    def m_dot_dyer(self, P_cc: float):
        """Return Dyer model mass flow rate.

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
