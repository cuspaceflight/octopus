from octopus import Fluid, thermo, FluidMixture, Manifold


def test_octopus_Fluid():
    water = Fluid("water")
    assert isinstance(water.chem, thermo.chemical.Chemical)


def test_octopus_FluidMixture():
    IPA95 = FluidMixture("isopropanol", "water", 0.95, 0.05)
    assert isinstance(IPA95.chem, thermo.chemical.Mixture)


def test_octopus_Manifold():
    IPA95 = FluidMixture("isopropanol", "water", 0.95, 0.05)
    manifold = Manifold(IPA95, 300, 1.75E6)
    assert isinstance(manifold.fluid.chem, thermo.chemical.Mixture)
