import thermo

class Fluid:
    """The fluid class is used to implement any fluids passing through the injector."""

    def __init__(self, name):
        self.name = name
        self.chem = thermo.chemical.Chemical(self.name)

class FluidMixture:
    """Easily create mixtures from Fluid objects, e.g. to observe the performance impact of additives."""

class Manifold:
    """The manifold class is used to organise injection elements by their propellant, and provide fluid property attributes to elements."""
