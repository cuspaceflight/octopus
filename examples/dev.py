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

    T = 250
    rhol_s = nitrous.rho_l(T)
    rhog_s = nitrous.rho_g(T)
    rho_range = linspace(rhog_s * 0.95, rhol_s * 1.05, 100)
    pv = [nitrous.get_properties(r, T)['p'] for r in rho_range]
    pp = [nitrous.p(r, T) for r in rho_range]
    plt.plot(rho_range, pv, color='green')
    plt.plot(rho_range, pp, color='red')
    plt.scatter(nitrous.rho_l(T), nitrous.p(nitrous.rho_l(T), T))
    plt.scatter(nitrous.rho_l(T), nitrous.VaporPressure(T), )
    plt.show()


if __name__ == "__main__":
    main()
