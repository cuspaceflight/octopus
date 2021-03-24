"""Implementation of injector classes.

References:
    - [1] - Thermophyiscal properties of nitrous oxide,
            IHS ESDU, http://edge.rit.edu/edge/P07106/public/Nox.pdf
"""

import thermo


class Fluid:
    """The fluid class is used to store fluid properties in the injector.
       This class is pretty useless at the moment, as thermo is the only
       option for properties - it just standardises behaviour between
       mixtures and pure fluids. If more chemical modules / sources are added,
       attributes will be standardised here.

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

    def __init__(self, name_1, source="thermo", **kwargs):
        self.source = source

        if self.source != "thermo":
            raise ValueError("Only \"thermo\" is supported currently.")

        if len(kwargs) % 2 != 1 and len(kwargs) > 0:
            raise ValueError("""Invalid number of kwargs: A name and mass
                             fraction required for each chemical""")

        self.name_list = []
        self.mf_list = []

        for kwarg in kwargs:
            if "name_" in kwarg:
                self.name_list.append(kwargs.get(kwarg))
                self.mf_list.append(kwargs.get(f"mf_{kwarg[5]}"))

        self.name_list.insert(0, name_1)
        self.mf_list.insert(0, kwargs.get("mf_1"))

        if len(self.name_list) != 1:
            if sum(self.mf_list) != 1:
                raise ValueError("Mass fractions provided do not sum to 1.")
        else:
            self.mf_list = [1]
            # Set mass fraction for a single element

        self.fluid = thermo.chemical.Mixture(self.name_list, ws=self.mf_list)
        # Configure the fluid object, which is a
        # mixture type even for only 1 chemical

    def rho(self, T, p):
        """Calculates fluid density under given conditions.
           Weighted average using mass fractions for mixtures.

        Args:
            T (float): Fluid temperature, Kelvin
            p (float): Fluid pressure, Pascals

        Returns:
            float: Fluid density, kg/m^3
        """

        rho = 0
        for i in range(len(self.name_list)):
            rho += self.mf_list[i] * thermo.chemical.Chemical(
                self.name_list[i], T=T, P=p).rho

        return rho

    def Cp(self, T, p):
        """Calculates specific heat capacity under given conditions.
           Weighted average using mass fractions for mixtures.

        Args:
            T (float): Fluid temperature, Kelvin
            p (float): Fluid pressure, Pascals

        Returns:
            float: Specific heat capacity, J/kg/K
        """

        Cp = 0
        for i in range(len(self.name_list)):
            Cp += self.mf_list[i] * thermo.chemical.Chemical(
                self.name_list[i], T=T, P=p).Cp

        return Cp

    def h(self, T, p):  # Ref [1], used for liquid phase enthalpy of nitrous
        pass

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
