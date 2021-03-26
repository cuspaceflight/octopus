"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""

import thermo


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
                orig_T = super().__getattribute__('T')  # temporarily saves T (T is an attribute, not a property)
                super().__setattr__('T', 10)  # temporarily sets T=10 to allow solver to reset
                super().Cpl  # computes Cpl to reset thermo's solver
                super().__setattr__('T', orig_T)  # resets T
                return super().Cpl  # returns value of Cpl using correct solver
        else:
            return temp_Cpl

    def h(self, T, p):  # Ref [1], used for liquid phase enthalpy of nitrous
        """DEBUG ONLY - Returns the value of h of nitrous oxide at -25C"""
        return -355000.0  # IMPLEMENT

    # Not sure what the best way to do this is - would be better if we
    # didn't have to manually reference tables on a per compound basis
    # thermopy3 is a candidate to do this with a different module

    def liqCheck(self, T, p):
        """Checks if the fluid contains any non-liquids.

        Args:
            T (float): Fluid temperature, Kelvin
            p (float): Fluid pressure, Pascals

        Returns:
            bool: True if all fluid components are liquid, else false.
        """

        for name in self.name_list:
            if thermo.chemical.Chemical(name, T=T, P=p).phase != "l":
                return False

        return True


class Manifold:
    """The manifold class is used to organise injection elements by
       their propellant, and provide fluid property attributes to elements."""

    def __init__(self, fluid, T_inlet, p_inlet):
        self.T_inlet = T_inlet  # Stagnation temperature, Kelvin
        self.p_inlet = p_inlet  # Stagnation, Pa
        self.fluid = fluid

        if self.fluid.liqCheck(T=T_inlet, p=p_inlet) is False:
            raise ValueError("Fluid is entering the manifold as a non-liquid!")

        self.rho = self.fluid.rho(T=T_inlet, p=p_inlet)
        self.Cp = self.fluid.Cp(T=T_inlet, p=p_inlet)


class Orifice:  # WIP
    """The orifice class is used to model thermodynamic changes in the
    fluid as it moves from the manifold into the combustion chamber """

    def __init__(self, fluid: thermo.chemical.Chemical, orifice_type, L, D):
        # subscript o represents initial conditions at stagnation
        self.fluid = fluid
        self.T_o = fluid.T
        self.p_o = fluid.P
        self.C_fo = fluid.Cpl
        self.v_fo = 1 / fluid.rhol
        self.v_fgo = abs(1 / fluid.rhol - 1 / fluid.rhog)
        self.h_fgo = abs(fluid.Hvap)
        self.orifice_type = orifice_type

        if orifice_type == 0:
            # omega is a parameter from Juang's paper relating pressure and temperature in the isentropic expansion
            self.omega = self.C_fo * self.T_o * self.p_o * (self.v_fgo / self.h_fgo) ** 2 / self.v_fo

    def v(self, p):
        if self.orifice_type == 0:
            # calculate specific volume from pressure using omega
            return self.v_fo * (self.omega * (self.p_o / p - 1) + 1)
