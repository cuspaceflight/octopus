import matplotlib.pyplot as plt

from octopus import Fluid, Orifice, PropertySource, Manifold


def main():
    nitrous = Fluid('nitrous oxide', P=18e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, T=400, method='thermo')
    nitrous_manifold = Manifold(nitrous, PropertySource(p=18e5, T=250))
    ipa_manifold = Manifold(isopropanol, PropertySource(p=18e5, T=400))
    nitrous_orifice = Orifice(nitrous_manifold, 1e-2, 1e-3)
    ipa_orifice = Orifice(ipa_manifold, 1e-2, 1e-3)

    plt.show()


if __name__ == "__main__":
    main()
