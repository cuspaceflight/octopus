import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from octopus import Fluid, Orifice, PropertySource, Manifold, Element


def main():
    nitrous = Fluid('nitrous oxide', P=20e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, T=400, method='thermo')
    nitrous_manifold = Manifold(nitrous, PropertySource(p=18e5, T=250))
    ipa_manifold = Manifold(isopropanol, PropertySource(p=18e5, T=400))
    nitrous_orifice = Orifice(nitrous_manifold, 1e-2, 2e-3)
    ipa_orifice = Orifice(ipa_manifold, 1e-2, 1e-3)
    element = Element([nitrous_orifice, nitrous_orifice], [ipa_orifice, ipa_orifice])

    print(element.of_ratio(16e5))

    T = 293
    p = 55e5

    x0 = [600, T]

    u = ['chi', 'p']
    y = [0, p]

    properties = least_squares(nitrous.fun_ps, x0, args=(u, y))
    rho, T = properties.x

    p = nitrous.get_properties(rho, T)['p']
    chi = nitrous.get_properties(rho, T)['chi']

    print(rho, p, T, chi)
    plt.show()


if __name__ == "__main__":
    main()
