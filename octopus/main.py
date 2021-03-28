"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""

from numpy import pi, sqrt, array, log
from scipy.optimize import least_squares
from thermo import chemical, property_molar_to_mass, utils


class Fluid(chemical.Chemical):
    """Inherits the thermo Chemical class, and represents a fluid with some overriden properties"""

    def __init__(self, ID: str, T: float = 298.15, P: float = 101325):
        super().__init__(ID, T, P)
        self.Cpl

    @property
    def Vl(self):
        return 1 / self.rhol

    @property
    def Vg(self):
        return 1 / self.rhog

    @property
    def Cpl(self):
        temp_Cpl = super().Cpl
        if self.ID == 'N2O':
            if temp_Cpl > 1500:  # If Cpl is below 1500, the supercritical polynomial is being used mistakenly by the solver
                return temp_Cpl
            else:
                orig_T = self.T  # temporarily saves T (T is an attribute, not a property)
                self.T = 10  # temporarily sets T=10 to allow solver to reset
                super().Cpl  # computes Cpl(10) to reset thermo's solver
                self.T = orig_T  # resets T
                return super().Cpl  # returns value of Cpl using correct solver
        else:
            return temp_Cpl

    @property
    def Hl(self):
        if self.ID == 'N2O':
            if not (self.Tm < self.T < self.Tc):
                raise Exception("Temperature must be within range 183-308K")
            Tr = self.T / self.Tc

            b1 = -200
            b2 = 116.043
            b3 = -917.225
            b4 = 794.779
            b5 = -589.587

            return 100 * (
                    b1 + b2 * (1 - Tr) ** (1 / 3) + b3 * (1 - Tr) ** (2 / 3) + b4 * (1 - Tr) + b5 * (1 - Tr) ** (4 / 3))
        else:
            raise Exception('Only liquid H data for N2O available')

    @property
    def Hg(self):
        if self.ID == 'N2O':
            if not (self.Tm < self.T < self.Tc):
                raise Exception("Temperature must be within range 183-308K")
            Tr = self.T / self.Tc
            b1 = -200
            b2 = 440.055
            b3 = -459.701
            b4 = 434.081
            b5 = -485.338

            return 100 * (
                    b1 + b2 * (1 - Tr) ** (1 / 3) + b3 * (1 - Tr) ** (2 / 3) + b4 * (1 - Tr) + b5 * (1 - Tr) ** (4 / 3))
        else:
            raise Exception('Only gaseous H data for N2O available')

    def calc_chi(self, M, V):
        v = V / M
        return (v - self.Vl) / (self.Vg - self.Vl)


class Manifold:
    """The manifold class is used to organise injection elements by
       their propellant, and provide fluid property attributes to elements."""

    def __init__(self, fluid: Fluid, T_inlet, p_inlet):
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

    def __init__(self, fluid: Fluid, chi0, P_cc, orifice_type, L, D, Cd=0.7):
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
        T, v = x
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
