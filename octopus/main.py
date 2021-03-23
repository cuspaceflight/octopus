import thermo


class Fluid:
    """The fluid class is used to store fluid properties in the injector.
       This class is pretty useless at the moment, as thermo is the only
       option for properties. If more modules / sources are added,
       methods will be implemented to standardise attributes.

        Attributes:
            name (string): Name of a chemical recognised by the source.
            source (string): Property source. Currently, only option is thermo.
                             Defaults to "thermo".
        """

    def __init__(self, name, source="thermo"):
        self.name = name

        if source == "thermo":
            self.chem = thermo.chemical.Chemical(self.name)
        else:
            raise ValueError("Only \"thermo\" is supported currently.")


class FluidMixture:
    """Easily create mixtures from Fluid objects, e.g. to observe
       the performance impact of additives.
       Currently just a container for the thermo equivalent, like Fluid.

       Might add option for any number of components with **kwargs
       Or a mixture of mixtures?

       Attributes:
            component_1 (string): Name of a chemical component of the
                                  mixture, recognised by the source.
            component_2 (string): Name of a chemical component of the
                                  mixture, recognised by the source.
            mf_1 (float): Mass fraction of component_1 in mixture.
            mf_2 (float): Mass fraction of component_2 in mixture.
            source (string): Property source. Currently, only option is thermo.
                             Defaults to "thermo".
    """

    def __init__(self, component_1, component_2, mf_1,
                 mf_2, source="thermo"):

        self.names = [component_1, component_2]

        if source == "thermo":
            self.chem = thermo.chemical.Mixture(
                                               [component_1, component_2],
                                               ws=[mf_1, mf_2])
        else:
            raise ValueError("Only \"thermo\" is supported currently.")


class Manifold:
    """The manifold class is used to organise injection elements by
       their propellant, and provide fluid property attributes to elements."""

    def __init__(self, fluid, T_inlet, p_inlet):

        self.T_inlet = T_inlet  # Kelvin
        self.p_inlet = p_inlet  # Pa
        self.fluid = fluid

        if isinstance(fluid, Fluid):
            self.fluid.chem.calculate(T=self.T_inlet, P=self.p_inlet)
            # Update thermo conditions

            if self.fluid.chem.phase != "l":
                print(self.fluid.chem.phase)
                raise Exception("""Fluid \"{}\" is entering the manifold
                                as a non-liquid!""".format(self.fluid.name))
            # Phase check at inlet

            self.rho = self.fluid.chem.rho
            # self.transport
            # Might use transport properties for manifold modelling later

        elif isinstance(fluid, FluidMixture):
            print("It's a mixture!")
            self.fluid.chem.calculate(T=self.T_inlet, P=self.p_inlet)
            for name in self.fluid.names:
                component = Fluid(name)
                component.chem.calculate(T=self.T_inlet, P=self.p_inlet)
                if component.chem.phase != "l":
                    raise Exception("""Fluid \"{}\" is entering the manifold
                                    as a non-liquid!""".format(name))
            # Phase check for each mixture component, at inlet
