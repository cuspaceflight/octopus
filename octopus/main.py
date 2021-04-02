"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""

from numpy import pi, sqrt, array, log, exp
from scipy.optimize import least_squares
from thermo import chemical, property_molar_to_mass, utils

from .utils import derivative


class Fluid(chemical.Chemical):
    """Inherits the thermo Chemical class, and represents a fluid with some overriden properties"""

    def __init__(self, ID: str, T: float = 298.15, P: float = 101325):
        super().__init__(ID, T, P)

    def alpha_0(self, delta, tau):
        """calculates the reduced helmholz idea gas energy if the chemical is nitrus oxide"""
        # symbol|superscript|subscipt --> c_p^0 = cp0
        if self.ID == 'N2O':

            a1 = -4.4262736272
            a2 = 4.3120475243
            c0 = 3.5
            c1 = 0  # not needed
            c2 = 1  # not needed

            v = array([2.1769, 1.6145, 0.48393])
            u = array([879.0, 2372.0, 5447.0])

            alpha_0 = a1 \
                      + a2 * tau \
                      + log(delta) \
                      + (c0 - 1) * log(tau) \
                      - (c1 * self.Tc ** c2) / (c2 * (c2 + 1)) * tau ** (-c2) \
                      + sum(v * log(1 - exp(-u * tau / self.Tc)))
            return alpha_0

        else:
            return None

    def alpha_r(self, delta, tau):
        """calculates the reduced helmholz residual energy if the chemical is nitrus oxide"""
        # symbol|superscript|subscipt --> c_p^0 = cp0
        if self.ID == 'N2O':

            n1 = 0.88045
            n2 = -2.4235
            n3 = 0.38237
            n4 = 0.068917
            n5 = 0.00020367
            n6 = 0.13122
            n7 = 0.46032
            n8 = -0.0036985
            n9 = -0.23263
            n10 = -0.00042859
            n11 = -0.042810
            n12 = -0.023038

            alpha_r = n1 * delta * tau ** 0.25 + \
                      n2 * delta * tau ** 1.25 + \
                      n3 * delta * tau ** 1.5 + \
                      n4 * delta ** 3 * tau ** 0.25 + \
                      n5 * delta ** 7 * tau ** 0.875 + \
                      n6 * delta * tau ** 2.375 * exp(-delta) + \
                      n7 * delta ** 2 * tau ** 2 * exp(-delta) + \
                      n8 * delta ** 5 * tau ** 2.125 * exp(-delta) + \
                      n9 * delta * tau * 3.5 * exp(-delta ** 2) + \
                      n10 * delta * tau ** 6.5 * exp(-delta ** 2) + \
                      n11 * delta ** 4 * tau ** 4.75 * exp(-delta ** 2) + \
                      n12 * delta ** 2 * tau ** 12.5 * exp(-delta ** 3)
            return alpha_r
        else:

            return None

    def rho_g(self, T):
        A = [-878.637631, 357.7403382, 574.66924864, -972.57488693, 630.47903568]
        tau = 1 - T / self.Tc
        return self.rhoc + A[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(A))

    def rho_l(self, T):
        A = [899.61701036, 179.53626729, 857.66247459, -2160.66039529, 2029.1923931]
        tau = 1 - T / self.Tc
        return self.rhoc + A[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(A))

    def a0_t(self, delta, tau):
        return derivative(self.alpha_0, 1, delta, tau)

    def a0_tt(self, delta, tau):
        return derivative(self.a0_t, 1, delta, tau)

    def ar_d(self, delta, tau):
        return derivative(self.alpha_r, 0, delta, tau)

    def ar_t(self, delta, tau):
        return derivative(self.alpha_r, 1, delta, tau)

    def ar_dd(self, delta, tau):
        return derivative(self.ar_d, 0, delta, tau)

    def ar_dt(self, delta, tau):
        return derivative(self.ar_d, 1, delta, tau)

    def ar_tt(self, delta, tau):
        return derivative(self.ar_t, 1, delta, tau)

    def dp_drho(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific*T*(1+2*delta*self.ar_d(delta,tau)+delta*delta*self.ar_dd(delta,tau))

    def dg_drho(self,rho,T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific*T*(2*self.ar_d(delta,tau)+delta*self.ar_dd(delta,tau))

    def p(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return rho * self.R_specific * T * \
               (1 + delta * self.ar_d(delta, tau))

    def u(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau))

    def h(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               (1 + tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) + delta * self.ar_d(delta, tau))

    def s(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * \
               (tau * (self.a0_t(delta, tau) + self.ar_t(delta, tau)) - self.alpha_0(delta, tau) - self.alpha_r(delta, tau))

    def g(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * \
               (1 + self.alpha_0(delta, tau) + self.alpha_r(delta, tau) + delta * self.ar_d(delta, tau))

    def p_g_error(self, x, Ts):
        rho_g, rho_l = x
        return array([self.p(rho_g, Ts) - self.p(rho_l, Ts),
                      self.g(rho_g, Ts) - self.g(rho_l, Ts)])

    def p_g_sep_args(self, rho_g, rho_l, Ts):
        return array([self.p(rho_g, Ts) - self.p(rho_l, Ts),
                      self.g(rho_g, Ts) - self.g(rho_l, Ts)])

    def jac(self, x, Ts):
        res = array([[self.dp_drho(x[0],Ts),self.dp_drho(x[1],Ts)],
                     [self.dg_drho(x[0],Ts),self.dg_drho(x[1],Ts)]])

        return res

    def rho_sat(self, T):

        x = array([self.rho_g(T),self.rho_l(T)])

        res = least_squares(self.p_g_error, x, jac=self.jac, args=[T])
        print(res.fun)
        return res.x

    def get_properties(self, rho, T):
        rhog = self.rho_g(T)
        rhol = self.rho_l(T)
        p = self.p(rhog, T)
        chi = (self.rhog - rho) * ((self.rhol - rho) / (self.rhol - self.rhog))
        h = self.h(rhog, T) * chi + self.h(rhol, T) * (1 - chi)
        s = self.s(rhog, T) * chi + self.s(rhol, T) * (1 - chi)

    class Manifold:
        """The manifold class is used to organise injection elements by
           their propellant, and provide fluid property attributes to elements."""

        def __init__(self, fluid, T_inlet, p_inlet):
            self.T_inlet = T_inlet  # Stagnation temperature, Kelvin
            self.p_inlet = p_inlet  # Stagnation, Pa
            self.fluid = fluid

            if self.fluid.phase != 'l':
                raise ValueError("Fluid is entering the manifold as a non-liquid!")

            self.rho = self.fluid.rho
            self.Cp = self.fluid.Cp

    class Orifice:  # WIP
        """The orifice class is used to model thermodynamic changes in the
        fluid as it moves from the manifold into the combustion chamber """

        def __init__(self, fluid, chi0, P_cc, orifice_type, L, D, Cd=0.7):
            # subscript o represents initial conditions at stagnation
            self.fluid = fluid
            self.P_o = fluid.P
            self.chi0 = chi0

            self.P_cc = P_cc
            self.orifice_type = orifice_type
            self.L = L
            self.D = D
            self.A = 0.2 * pi * self.D ** 2
            self.Cd = Cd

            self.c_lo = fluid.Cpl
            self.h_o = chi0 * fluid.Hg + (chi0 - 1) * fluid.Hl
            self.s_o = None

            self.v_lgo = abs(1 / fluid.rhol - 1 / fluid.rhog)
            self.h_lgo = abs(fluid.Hvap)

            self.T_o = fluid.Tsat(self.P_o)
            self.fluid.T = self.T_o
            self.v_o = fluid.Vl + chi0 * self.v_lgo

            self.h_1 = 0

            '''
            if orifice_type == 0:
                # omega is a parameter from Juang's paper relating pressure and temperature in the isentropic expansion
                # eta is the critical pressure ratio
                self.omega = self.C_fo * self.T_o * self.P_o * (self.v_fgo / self.h_fgo) ** 2 / self.v_fo
                self.eta = 0.55 + 0.217 * log(self.omega) - 0.046 * log(self.omega) ** 2 + 0.004 * log(self.omega) ** 3
            '''

        def m_dot_SPI(self):
            if self.orifice_type == 0:
                return self.Cd * self.A * sqrt(2 * 1 / self.v_o * (self.P_o - self.P_cc))

        def f(self, x, u):

            chi, fluid = u
            P = fluid.VaporPressure(T)

            # liquid entropy at outlet wrt critical pressure (equal for liquid & gas)
            Slm = fluid.HeatCapacityLiquid.T_dependent_property_integral_over_T(309, T) + fluid.eos.to_TP(T,
                                                                                                          P).S_dep_l - utils.R * log(
                P / fluid.Pc)
            Sl = property_molar_to_mass(Slm, fluid.MW)
            # gaseous entropy at outlet wrt critical pressure (equal for liquid & gas)
            Sgm = fluid.HeatCapacityGas.T_dependent_property_integral_over_T(309, T) + fluid.eos.to_TP(T,
                                                                                                       P).S_dep_g - utils.R * log(
                P / fluid.Pc)
            Sg = property_molar_to_mass(Sgm, fluid.MW)

            Vl = property_molar_to_mass(fluid.eos.to_TP(T, P).V_l, fluid.MW)
            Vg = property_molar_to_mass(fluid.eos.to_TP(T, P).V_g, fluid.MW)

            p_res = P
            chi_res = (v - Vl) / (Vg - Vl)
            S_res = Sl + chi_res * (Sg - Sl)

            return array([p_res, S_res])

        def fun(self, x, u, y):
            return self.f(x, u) - y

        def m_dot_HEM(self):
            if self.orifice_type == 0:
                P_o = self.P_o
                T_o = self.fluid.Tsat(P_o)
                # ambient pressure at outlet = P_cc
                P_1 = self.P_cc
                # liquid entropy at inlet
                S_1m = self.fluid.HeatCapacityLiquid.T_dependent_property_integral_over_T(self.fluid.Tc, T_o) \
                       + self.fluid.eos.to_TP(T_o, P_o).S_dep_l - utils.R * log(P_o / self.fluid.Pc)
                S_1 = property_molar_to_mass(S_1m, self.fluid.MW)

                print(f'set:{S_1}')

                # measured values are exit pressure and entropy (s=const)
                y = array([P_1, S_1])
                # independent variables are the fluid initial state
                u = array([self.chi0, self.fluid])
                # initial guess of T and v
                x0 = array([245, self.fluid.Vg / 2])
                # bounds on T
                bounds = ([184, self.fluid.Vl], [308, self.fluid.Vg])

                res = least_squares(self.fun, x0, args=(u, y), bounds=bounds)
                T_1, v_1 = res.x
                P, S = self.f([T_1, v_1], [0, self.fluid])

                # set fluid
                self.fluid.calculate(T_1, P_1)

                Vl = property_molar_to_mass(self.fluid.eos.to_TP(T_1, P_1).V_l, self.fluid.MW)
                Vg = property_molar_to_mass(self.fluid.eos.to_TP(T_1, P_1).V_g, self.fluid.MW)
                chi_1 = (v_1 - Vl) / (Vg - Vl)

                print(f'chi: {chi_1}')
                print(f'T: {bounds[0][0]},\t\t\t\t\t\t {T_1}, \t{bounds[1][0]}')
                print(f'v: {bounds[0][1]}, \t{v_1}, \t{bounds[1][1]}')

                print(P, S)

                # calculate dh
                self.h_1 = chi_1 * self.fluid.Hg + (1 - chi_1) * self.fluid.Hl

                # print(chi_1, v_1, self.fluid.Vg, self.fluid.eos.to_TP(T_1, self.P_cc).V_g)

                return self.Cd * self.A * sqrt(self.h_o - self.h_1) / v_1

        def m_dot_dyer(self):
            kappa = 1  # fluid enters orifice at vapour pressure
            return 1 / (1 + kappa) * self.m_dot_HEM() + (1 - 1 / (1 + kappa)) * self.m_dot_SPI()
