from thermo import chemical

from octopus import Fluid, Manifold


def test_octopus_Fluid():
    water = Fluid("water")
    assert isinstance(water, chemical.Chemical)


'''
def test_octopus_FluidMixture():
    IPA95 = Fluid("isopropanol", name_2="water", mf_1=0.95, mf_2=0.05)
    assert isinstance(IPA95.fluid, thermo.chemical.Mixture)
'''


def test_octopus_Manifold():
    IPA = Fluid("isopropanol")
    manifold = Manifold(IPA, 300, 1.75E6)
    assert isinstance(manifold.fluid, chemical.Chemical)
