"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""

import thermo
from numpy import pi, sqrt, log


class Fluid(thermo.chemical.Chemical):
    """Inherits the thermo Chemical class, and represents a fluid with some overriden properties

       Keyword Args:
            name_n (string): Name of the nth mixture component.
            mf_n (float): Mass fraction of nth component.
                          The sum of all mass fractions must equal 1.

        Attributes:
            fluid (type depends on source): Mixture or pure fluid object.
            source (string): Property source. Currently, only option is thermo.
                             Defaults to "thermo".
            name_list (list): List of chemical names, same order as mf_list.
            mf_list (list): List of mass fractions, same order as names_list.
        """

    def __init__(self, ID: str, T: float = 298.15, P: float = 101325):
        super().__init__(ID, T, P)

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

    def __init__(self, fluid: thermo.chemical.Chemical, orifice_type, L, D):
        # subscript o represents initial conditions at stagnation
        self.fluid = fluid
        self.T_o = fluid.T
        self.P_o = fluid.P
        self.C_fo = fluid.Cpl
        self.v_fo = 1 / fluid.rhol
        self.v_fgo = abs(1 / fluid.rhol - 1 / fluid.rhog)
        self.h_fgo = abs(fluid.Hvap)
        self.orifice_type = orifice_type
        self.L = L
        self.D = D

        if orifice_type == 0:
            # omega is a parameter from Juang's paper relating pressure and temperature in the isentropic expansion
            # eta is the critical pressure ratio
            self.omega = self.C_fo * self.T_o * self.P_o * (self.v_fgo / self.h_fgo) ** 2 / self.v_fo
            self.eta = 0.55 + 0.217 * log(self.omega) - 0.046 * log(self.omega) ** 2 + 0.004 * log(self.omega) ** 3

    def v(self, P):
        if self.orifice_type == 0:
            # calculate specific volume from pressure using omega
            return self.v_fo * (self.omega * (self.P_o / P - 1) + 1)

    def m_dot(self):
        if self.orifice_type == 0:
            # m_dot = A*sqrt(Po/vo)*eta/sqrt(omega)
            return (pi * self.D ** 2) * self.eta \
                   * sqrt(self.P_o / self.v_fo / self.omega)
