import matplotlib.pyplot as plt
from numpy import linspace

from octopus import Fluid, Orifice, PropertySource, Manifold, Element


def main():
    nitrous = Fluid('nitrous oxide', P=20e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, T=400, method='thermo')
    nitrous_manifold = Manifold(nitrous, PropertySource(p=18e5, T=250))
    ipa_manifold = Manifold(isopropanol, PropertySource(p=18e5, T=400))
    nitrous_orifice = Orifice(nitrous_manifold, 1e-2, 2e-3)
    ipa_orifice = Orifice(ipa_manifold, 1e-2, 1e-3)
    element = Element([nitrous_orifice, nitrous_orifice], [ipa_orifice, ipa_orifice])

    pp = linspace(16.7e5 - 1e5, 16.7e5 + 1e5, 40)
    of = [element.of_ratio(p) for p in pp]

    plt.plot(pp, of, label='O/F')
    plt.xlabel('chamber pressure (Pa)')
    plt.show()


if __name__ == "__main__":
    main()
