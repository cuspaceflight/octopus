"""Implementation of injector fluid dynamics classes.

"""
from dataclasses import dataclass
from typing import Iterable, Sequence

import CoolProp.CoolProp as CP
import numpy as np
import scipy.optimize
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
        x = self.chi()
        return 1 / (1 / (1 + S * ((1 - x) / x) * self.rhog([self.state.T()])[0] / self.rhol([self.state.T()])[0]))

    def chi(self):
        x = self.state.Q()
        if x >= 0:
            return x
        if self.state.rhomass() > self.rhol([self.state.T()])[0]:
            return 0
        else:
            return 1

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
        T = self.state.T()
        mug = self.viscosityg(T)
        mul = self.viscosityl(T)
        x = self.chi()
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
    pos: float
    p: float
    T: float
    v: float
    rho: float
    h: float
    mu: float
    A: float
    D: float
    s: float

    def __repr__(self):
        return f'Cell[p: {self.p:.0f}, v: {self.v:2.1f}, rho: {self.rho:.3f}, mu: {self.mu:.3f}, A: {self.A:.5f}]'


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

    def __init__(self, fluid: Fluid, parent: PropertySource, A: float):
        """Initialise :class:`Manifold` object with a working fluid and property parent.

        :param fluid: :class:`Fluid` object to use EOS functions from
        :param parent: :class:`PropertySource` object to get p and T from
        """
        self.fluid = fluid
        self.parent = parent
        self.chi = 0
        self.A = A

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
    ANNULAR = 2

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
                print('no pressure drop in SPI calc')
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
                print('no pressure drop in HEM calc')

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
                print('no decrease in enthalpy in HEM calc')

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

    def sub(self, pcc, mdot):
        return mdot - self.m_dot_dyer(pcc)

    def p_dyer(self, mdot):

        res = scipy.optimize.least_squares(self.sub, 13e5, args=[mdot], bounds=[1e5, 20e5], ftol=0.0001, gtol=None)
        # print(mdot, self.m_dot_dyer(res.x)[0])
        return res.x[0]

    def p_patel_dyer(self, mdot):
        # aliases
        fluid = self.manifold.fluid
        fluid.set_state(P=self.manifold.p, T=self.manifold.T)
        # run dyer to get inlet conditions
        s0 = fluid.state.smass()
        p0 = self.p_dyer(mdot)
        print(f'inlet p: {p0 / 1e5:.1f} bar\n'
              f'inlet s: {s0:.0f} kJ/kg.K')
        fluid.set_state(P=p0, S=s0)

        # setup distributions
        N = 100
        pos_dist = np.linspace(0, self.L, N)
        A = self.A * np.ones_like(pos_dist)
        Dh = self.D * np.ones_like(pos_dist)
        p = np.interp(pos_dist, [0, self.L], [p0, 0.9 * p0])

        # init cells
        cells = []
        for i in range(N):
            cells.append(Cell(pos=pos_dist[i], p=p[i], A=A[i], D=Dh[i], T=0, v=0, rho=0, h=0, s=None, mu=fluid.viscosity()))

        # known values in first cell
        cells[0].T = self.manifold.T
        cells[0].mu = fluid.viscosity()
        cells[0].rho = fluid.state.rhomass()
        cells[0].h = fluid.state.hmass()
        cells[0].v = mdot / (cells[0].rho * cells[0].A)

        # run advance repeatedly on cells
        sets_0 = [[cell.p for cell in cells]]
        sets_1 = [[cell.v for cell in cells]]
        sets_2 = [[cell.s for cell in cells]]
        sets_3 = [[cell.rho for cell in cells]]
        for i in range(10):
            print(i)
            self.p_patel_advance(cells, mdot)
            sets_0.append([cell.p for cell in cells])
            sets_1.append([cell.v for cell in cells])
            sets_2.append([cell.s for cell in cells])
            sets_3.append([cell.rho for cell in cells])

        return sets_0, sets_1, sets_2, sets_3

    def p_patel(self, mdot):
        # aliases
        fluid = self.manifold.fluid
        fluid.set_state(P=self.manifold.p, T=self.manifold.T)

        # setup distributions
        N = 1000
        pos_dist = np.linspace(0, 2 * self.L, N)
        A = np.interp(pos_dist, [0, self.L, 2 * self.L], [self.manifold.A, self.A, self.A])
        inlet_Dh = np.sqrt(4 * self.manifold.A / np.pi)
        Dh = np.interp(pos_dist, [0, self.L, 2 * self.L], [inlet_Dh, self.D, self.D])
        p = np.interp(pos_dist, [0, self.L], [self.manifold.p, 0.9 * self.manifold.p])

        # init cells
        cells = []
        for i in range(N):
            cells.append(Cell(pos=pos_dist[i], p=p[i], A=A[i], D=Dh[i], T=0, v=0, rho=0, h=0, s=None, mu=fluid.viscosity()))

        # known values in first cell
        cells[0].T = self.manifold.T
        cells[0].mu = fluid.viscosity()
        cells[0].rho = fluid.state.rhomass()
        cells[0].h = fluid.state.hmass()
        cells[0].v = mdot / (cells[0].rho * cells[0].A)

        # run advance repeatedly on cells
        sets_0 = [[cell.p for cell in cells]]
        sets_1 = [[cell.v for cell in cells]]
        sets_2 = [[cell.s for cell in cells]]
        sets_3 = [[cell.rho for cell in cells]]
        for i in range(5):
            print(i)
            self.p_patel_advance(cells, mdot)
            sets_0.append([cell.p for cell in cells])
            sets_1.append([cell.v for cell in cells])
            sets_2.append([cell.s for cell in cells])
            sets_3.append([cell.rho for cell in cells])

        return sets_0, sets_1, sets_2, sets_3

    def p_patel_advance(self, cells, mdot):
        fluid = self.manifold.fluid
        # initialise fluid
        fluid.set_state(P=cells[0].p, T=cells[0].T)
        # calculate derivative of pressure with respect to position
        dx_l = np.diff([cell.pos for cell in cells])
        dp_dx_l = np.diff([cell.p for cell in cells]) / dx_l
        # calculate h+0.5*v**2 from first cell
        h_const = cells[0].h + 0.5 * cells[0].v ** 2

        for cell in cells:
            try:
                # changes from last iteration
                i = cells.index(cell)
                dp_dx = dp_dx_l[min(i, len(dp_dx_l) - 1)]
                dx = dx_l[min(i, len(dp_dx_l) - 1)]

                # calculate density and h
                cell.rho = mdot / (cell.v * cell.A)
                cell.h = h_const - 0.5 * cell.v ** 2
                # calculate viscosity
                fluid.set_state(P=cell.p, H=cell.h)
                cell.mu = fluid.viscosity()

                # calculate kinematic properties
                Re = mdot * cell.D / (cell.mu * cell.A)

                # calculate v for next cell
                dv_dx = -(-fd(Re) * 0.5 * cell.rho * cell.v ** 2 * 4 / cell.D + dp_dx) / (cell.rho * cell.v)
                cells[min(i + 1, len(cells) - 1)].v = cell.v + dv_dx * dx

                # calculate new pressure from density and enthalpy
                fluid.set_state(D=cell.rho, H=cell.h)
                cell.p = fluid.state.p()
                try:
                    cell.s = fluid.state.smass()
                except ValueError:
                    cell.s = None
            except TypeError:
                pass

        return cells

    def p_collins(self, p_cc):
        # aliases
        fluid = self.manifold.fluid
        fluid.set_state(P=self.manifold.p, T=self.manifold.T)

        # setup distributions
        N = 1000
        pos_dist = np.linspace(0, 2 * self.L, N)
        A = np.interp(pos_dist, [0, 1 * self.L, 2 * self.L], [self.manifold.A, self.A, self.A])
        inlet_Dh = np.sqrt(4 * self.manifold.A / np.pi)
        Dh = np.interp(pos_dist, [0, 1 * self.L, 2 * self.L], [inlet_Dh, self.D, self.D])
        p = np.interp(pos_dist, [0, 2 * self.L], [self.manifold.p, 0.2 * self.manifold.p])

        # init cells
        cells = []
        for i in range(N):
            cells.append(Cell(pos=pos_dist[i], p=p[i], A=A[i], D=Dh[i], T=0, v=0, rho=0, h=0, s=0, mu=fluid.viscosity()))

            # known values in first cell
            cells[0].T = self.manifold.T
            cells[0].mu = fluid.viscosity()
            cells[0].rho = fluid.state.rhomass()
            cells[0].h = fluid.state.hmass()

    def dP_HEM_frictional(self, mdot: float):
        # p0,T0 known
        p0 = self.manifold.p
        T0 = self.manifold.T
        fluid = self.manifold.fluid

        # find initial conditions
        fluid.set_state(P=p0, T=T0)
        s0 = fluid.state.smass()

        # calculate inlet pressure drop
        p_inlet = self.p_dyer(mdot)

        # set conditions after inlet
        fluid.set_state(S=s0, P=p_inlet)

        # setup integration values
        num_cells = 500
        dz = self.L / num_cells

        # constant properties
        G = mdot / self.A

        # cell arrays for integration
        cells = [Cell(p=p_inlet, x=fluid.chi())] + [Cell(p=p_inlet, x=fluid.chi()) for _ in range(num_cells)]
        print(fluid.state.T())

        try:
            for i in range(num_cells - 1):
                vf = 1 / fluid.rhol([fluid.state.T()])[0]
                vfg = 1 / fluid.rhog([fluid.state.T()])[0] - vf

                fluid.set_state(P=cells[i].p, S=s0)
                cells[i + 1].x = fluid.chi()

                dx_dz = np.clip((cells[i + 1].x - cells[i].x) / dz, -0.3, 0.3)
                dvg_dP = fluid.dvg_dP_saturation()

                V = G / fluid.state.rhomass()
                Re = fluid.state.rhomass() * V * self.D / fluid.viscosity()

                dP_dz = -((2 * (fd(Re) / 4) * G ** 2 * vf / self.D) * (1 + cells[i + 1].x * (vfg / vf)) + G ** 2 * vf * (
                        vfg / vf) * dx_dz) / (1 + G ** 2 * cells[i + 1].x * dvg_dP)

                cells[i + 2].p = cells[i].p + dP_dz * dz
                print(f'{dP_dz * dz:7.0f}\t{dx_dz:4.0f}\t{V:3.0f}')

        except Exception as e:
            print('intgration failed in dP_HEM_frictional')
            print(e)

        return [cell.p / 100000 for cell in cells]


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
