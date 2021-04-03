"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
    - [2] - An investigation of injectors for use with high vapor pressure
            propellants with applications to hybrid rockets, Benjamin S. Waxman
            https://stacks.stanford.edu/file/druid:ng346xh6244/BenjaminWaxmanFinal-augmented.pdf
    - [3] - Short Fundamental Equations of State for 20 Industrial Fluids,
            Lemmon and Span, https://pubs.acs.org/doi/pdf/10.1021/je050186n
"""

from numpy import pi, sqrt, array, log, exp, float64, nan_to_num
from scipy.optimize import least_squares
from thermo import chemical


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
        return nan_to_num(self.rhoc + A[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(A)))

    def rho_l(self, T):
        A = [899.61701036, 179.53626729, 857.66247459, -2160.66039529, 2029.1923931]
        tau = 1 - T / self.Tc
        return nan_to_num(self.rhoc + A[0] * (tau ** 0.35 - 1) + sum(a * tau ** i for i, a in enumerate(A)))

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
        return self.R_specific * T * (1 + 2 * delta * self.ar_d(delta, tau) + delta * delta * self.ar_dd(delta, tau))

    def dg_drho(self, rho, T):
        delta = rho / self.rhoc
        tau = self.Tc / T
        return self.R_specific * T * (2 * self.ar_d(delta, tau) + delta * self.ar_dd(delta, tau))

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

    def cpl(self, T):
        """Specific heat at constant pressure for saturated liquid N2O.
           See Ref [1].

        Args:
            T (float): Nitrous oxide temperature, K

        Returns:
            (float): Cp for liquid N2O at given temperature.
        """
        if not 183.15 <= T <= 303.15:
            raise ValueError(f"Temperature ({T} K) out of range")
        Tr = 309.57  # Find the reduced temperature (T / T_crit)
        return 2.49973*(1 + 0.023454/(1-Tr) - 3.80136*(1-Tr) +
                        13.0945*(1-Tr)**2 - 14.5180*(1-Tr)**3)

    def z(self, delta, tau):
        return 1 + delta*ar_d(delta, tau)
        # Return the compressibility - see Ref [3], Eqn (15)

    # I might not be following your intended syntax here
    def dz_dT(self, ):
        return derivative(z, 0, )

    def fun_ps(self, x, u, y):
        return [self.get_properties(x[0], x[1])[var] - val for var, val in zip(u, y)]

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
            chi = float64(0.0)
            h = self.h(rho, T)
            s = self.s(rho, T)
        else:
            p = self.VaporPressure(T)
            chi = rhog * (rhol - rho) / (rho * (rhol - rhog))
            h = self.h(rhog, T) * chi + self.h(rhol, T) * (1 - chi)
            s = self.s(rhog, T) * chi + self.s(rhol, T) * (1 - chi)

        res = {'p': p, 'chi': chi, 'h': h, 's': s}
        for key in res:
            if isinstance(res[key], type(None)):
                print(key,res[key],rho,rhol,rhog)
        return res


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

        '''
        if orifice_type == 0:
            # omega is a parameter from Juang's paper relating pressure and temperature in the isentropic expansion
            # eta is the critical pressure ratio
            self.omega = self.C_fo * self.T_o * self.P_o * (self.v_fgo / self.h_fgo) ** 2 / self.v_fo
            self.eta = 0.55 + 0.217 * log(self.omega) - 0.046 * log(self.omega) ** 2 + 0.004 * log(self.omega) ** 3
        '''

    def m_dot_SPI(self):
        if self.orifice_type == 0:
            p0 = self.P_o
            chi0 = self.chi0
            u = ['p', 'chi']
            y = [p0, chi0]
            initial = least_squares(self.fluid.fun_ps, [800, 250], args=[u, y])
            rho0, T0 = initial.x
            return self.A * sqrt(2 * rho0 * (self.P_o - self.P_cc))

    def m_dot_HEM(self):
        if self.orifice_type == 0:
            # find initial conditions
            # chi0,p0 known
            p0 = self.P_o
            chi0 = self.chi0

            u = ['p', 'chi']
            y = [p0, chi0]
            initial = least_squares(self.fluid.fun_ps, [800, 250], args=[u, y])
            rho0, T0 = initial.x
            h0 = self.fluid.get_properties(rho0, T0)['h']
            s0 = self.fluid.get_properties(rho0, T0)['s']

            # set final conditions
            # p1, s1 known
            p1 = self.P_cc
            s1 = s0

            u = ['p', 's']
            y = [p1, s1]

            final = least_squares(self.fluid.fun_ps, [500, 245], args=[u, y])
            rho1, T1 = final.x
            h1 = self.fluid.get_properties(rho1, T1)['h']
            X1 = self.fluid.get_properties(rho1, T1)['chi']

            return self.A * rho1 * sqrt(h0 - h1)

    def m_dot_dyer(self):
        kappa = 1  # fluid enters orifice at vapour pressure
        W = 1 / (1 + kappa)
        return self.Cd * ((1 - W) * self.m_dot_SPI() + W * self.m_dot_HEM())

    def Y(self, Fluid):
        """Calculate the general compressibility correction factor for the orifice.
        See Ref [2], section 2.1.1.2.

        Returns:
            float: Mass flow compressibility correction.
        """
        T = Fluid.T
        P = Fluid.P
        R = 8.31446/(Fluid.MW/1000)  # Get the specific gas constant
        #  cpl = Fluid.cpl(T)  # Saturated liquid Cp at this temperature
        #  gamma = cpl/(cpl - R)  # Ratio of specific heats
        rho_l = Fluid.rho_l(T)

        # Isentropic power law exponent - equation 2.22 of Ref [2]
        Z = P/(rho_l*R*T)  # Compressibiliy
        #  dZdT_rho = None

        #  n = gamma * (Z + T)/(Z+T)

        return None
