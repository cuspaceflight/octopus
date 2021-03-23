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
       the performance impact of additives."""


class Manifold:
    """The manifold class is used to organise injection elements by
       their propellant, and provide fluid property attributes to elements."""
