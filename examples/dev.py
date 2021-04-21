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

    pl = []
    pg = []
    pv = []
    T=290
    rho = linspace(nitrous.rho_g(T)*0.5, nitrous.rho_l(T)*1.2, 100)
    for r in rho:
        pl.append(nitrous.p(r, T))
        pg.append(nitrous.p(r, T))
        pv.append(nitrous.VaporPressure(T))
    plt.plot(rho, pl, label='l')
    plt.plot(rho, pg, label='g')
    plt.plot(rho, pv, label='v')
    plt.scatter(nitrous.rho_g(T),nitrous.VaporPressure(T))
    plt.scatter(nitrous.rho_l(T), nitrous.VaporPressure(T))
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
