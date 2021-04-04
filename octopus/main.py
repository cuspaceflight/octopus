"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""
from functools import lru_cache
from json import load
from os.path import dirname

from numpy import pi, sqrt, array, log, exp, nan_to_num, inf
from scipy.optimize import least_squares
from thermo import chemical

from .utils import derivative

cwd = dirname(__file__)


class Fluid(chemical.Chemical):
    """Inherits the thermo Chemical class, and represents a fluid with some overriden properties"""

    def __init__(self, ID: str, T: float = 298.15, P: float = 101325):
        super().__init__(ID, T, P)

        with open(f'{cwd}\\{self.ID}.json', 'r') as f:
            data = load(f)

        self.a = data['a']
        self.c = data['c']
        self.v = array(data['v'])
        self.u = array(data['u'])
        self.n = data['n']
        self.Ag = data['Ag']
        self.Al = data['Al']
        self.polar = bool(data['polar'])

    @lru_cache(maxsize=1)
    def alpha_0(self, delta, tau):
        """Calculate the reduced ideal gas Helmholz energy is the fluid coefficients are known"""

        alpha_0 = (self.a[0]
                   + self.a[1] * tau
                   + log(delta)
                   + (self.c[0] - 1) * log(tau)
                   - (self.c[1] * self.Tc ** self.c[2]) / (self.c[2] * (self.c[2] + 1)) * tau ** (-self.c[2])
                   + sum(self.v * log(1 - exp(-self.u * tau / self.Tc))))
        return alpha_0

    @lru_cache(maxsize=1)
    def alpha_r(self, delta, tau):
        """calculates the reduced helmholz residual energy if the chemical is nitrus oxide"""

        if self.polar:

            alpha_r = (self.n[0] * delta * tau ** 0.25
                       + self.n[1] * delta * tau ** 1.25
                       + self.n[2] * delta * tau ** 1.5
                       + self.n[3] * delta ** 3 * tau ** 0.25
                       + self.n[4] * delta ** 7 * tau ** 0.875
                       + self.n[5] * delta * tau ** 2.375 * exp(-delta)
                       + self.n[6] * delta ** 2 * tau ** 2 * exp(-delta)
                       + self.n[7] * delta ** 5 * tau ** 2.125 * exp(-delta)
                       + self.n[8] * delta * tau * 3.5 * exp(-delta ** 2)
                       + self.n[9] * delta * tau ** 6.5 * exp(-delta ** 2)
                       + self.n[10] * delta ** 4 * tau ** 4.75 * exp(-delta ** 2)
                       + self.n[11] * delta ** 2 * tau ** 12.5 * exp(-delta ** 3))
            return alpha_r
        else:
            raise NotImplementedError

    def rho_g(self, T):
        tau = 1 - T / self.Tc
        return nan_to_num(self.rhoc + self.Ag[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(self.Ag)))

    def rho_l(self, T):
        tau = 1 - T / self.Tc
        return nan_to_num(self.rhoc + self.Al[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(self.Al)))

    @lru_cache(maxsize=1)
    def a0_t(self, delta, tau):
        return derivative(self.alpha_0, 1, delta, tau)

    @lru_cache(maxsize=1)
    def a0_tt(self, delta, tau):
        return derivative(self.a0_t, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_d(self, delta, tau):
        return derivative(self.alpha_r, 0, delta, tau)

    @lru_cache(maxsize=1)
    def ar_t(self, delta, tau):
        return derivative(self.alpha_r, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_dd(self, delta, tau):
        return derivative(self.ar_d, 0, delta, tau)

    @lru_cache(maxsize=1)
    def ar_dt(self, delta, tau):
        return derivative(self.ar_d, 1, delta, tau)

    @lru_cache(maxsize=1)
    def ar_tt(self, delta, tau):
        return derivative(self.ar_t, 1, delta, tau)

    @lru_cache(maxsize=1)
    def dp_drho(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * (1 + 2 * delta * self.ar_d(delta, tau) + delta * delta * self.ar_dd(delta, tau))

    @lru_cache(maxsize=1)
    def dg_drho(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * (2 * self.ar_d(delta, tau) + delta * self.ar_dd(delta, tau))

    @lru_cache(maxsize=1)
    def p(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return rho * self.R_specific * T * \
               (1 + delta * self.ar_d(delta, tau))

    @lru_cache(maxsize=1)
    def u(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau))

    @lru_cache(maxsize=1)
    def h(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               (1 + tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) + delta * self.ar_d(delta, tau))

    @lru_cache(maxsize=1)
    def s(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * \
               (tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) - self.alpha_0(delta, tau) - self.alpha_r(delta, tau))

    @lru_cache(maxsize=1)
    def g(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               (1 + self.alpha_0(delta, tau) + self.alpha_r(delta, tau) + delta * self.ar_d(delta, tau))

    def fun_ps(self, x, u, y):
        return [self.get_properties(x[0], x[1])[var] - val for var, val in zip(u, y)]

    @lru_cache(maxsize=1)
    def get_properties(self, rho, T):
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


class Manifold:
    """Organises Orifices into Elements, contains propellants"""

    def __init__(self, fluid, T_inlet, p_inlet):
        self.T_inlet = T_inlet  # Stagnation temperature, Kelvin
        self.p_inlet = p_inlet  # Stagnation, Pa
        self.fluid = fluid

        if self.fluid.phase != 'l':
            raise ValueError("Fluid is entering the manifold as a non-liquid!")

        self.rho = self.fluid.rho
        self.Cp = self.fluid.Cp


class Orifice:  # WIP
    """Model the thermodynamic changes as fluid moves through the orifice"""

    def __init__(self, fluid, chi0, orifice_type, L, D, Cd=0.7):
        # subscript o represents initial conditions at stagnation
        self.fluid = fluid
        self.P_o = fluid.P
        self.chi0 = chi0

        self.orifice_type = orifice_type
        self.L = L
        self.D = D
        self.A = 0.2 * pi * self.D ** 2
        self.Cd = Cd

    @lru_cache(maxsize=1)
    def m_dot_SPI(self, P_cc):
        if self.orifice_type == 0:
            # find initial conditions
            # chi0,p0 known
            p0 = self.P_o
            chi0 = self.chi0

            if p0 < P_cc: return None

            u = ['p', 'chi']
            y = [p0, chi0]

            initial = least_squares(self.fluid.fun_ps, [800, 250], bounds=([0, 0], [inf, self.fluid.Tc]), args=[u, y])
            rho0, T0 = initial.x

            return self.A * sqrt(2 * rho0 * (self.P_o - P_cc))

    @lru_cache(maxsize=1)
    def m_dot_HEM(self, P_cc):
        if self.orifice_type == 0:
            # find initial conditions
            # chi0,p0 known
            p0 = self.P_o
            chi0 = self.chi0

            if p0 < P_cc: return None

            u = ['p', 'chi']
            y = [p0, chi0]

            initial = least_squares(self.fluid.fun_ps, [800, 250], bounds=([0, 0], [inf, self.fluid.Tc]), args=[u, y])
            rho0, T0 = initial.x

            h0 = self.fluid.get_properties(rho0, T0)['h']
            s0 = self.fluid.get_properties(rho0, T0)['s']

            # set final conditions
            # p1, s1 known
            p1 = P_cc
            s1 = s0

            u = ['p', 's']
            y = [p1, s1]

            final = least_squares(self.fluid.fun_ps, [500, 245], bounds=([0, 0], [inf, self.fluid.Tc]), args=[u, y])
            rho1, T1 = final.x

            h1 = self.fluid.get_properties(rho1, T1)['h']
            X1 = self.fluid.get_properties(rho1, T1)['chi']
            if h1 > h0: return None
            return self.A * rho1 * sqrt(h0 - h1)

    @lru_cache(maxsize=1)
    def m_dot_dyer(self, P_cc):
        kappa = 1  # fluid enters orifice at vapour pressure
        W = 1 / (1 + kappa)
        if not (self.m_dot_SPI(P_cc) and self.m_dot_HEM(P_cc)): return None
        return self.Cd * ((1 - W) * self.m_dot_SPI(P_cc) + W * self.m_dot_HEM(P_cc))
