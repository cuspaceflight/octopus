import matplotlib.pyplot as plt
import numpy as np

from octopus import Fluid, Orifice, PropertySource, Manifold


def main():
    nitrous = Fluid('nitrous oxide', P=18e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, method='thermo')
    nitrous_manifold = Manifold(nitrous, PropertySource(p=18e5, T=250))
    ipa_manifold = Manifold(isopropanol, PropertySource(p=18e5, T=400))
    nitrous_orifice = Orifice(nitrous_manifold, 1e-2, 1e-3)
    ipa_orifice = Orifice(ipa_manifold, 1e-2, 1e-3)

    P_cc = np.linspace(0, 19e5, 50)
    SPI = []
    SPIn = []
    HEM = []
    DYER = []
    for Pcc in P_cc:
        SPI.append(ipa_orifice.m_dot_SPI(Pcc))
        SPIn.append(nitrous_orifice.m_dot_SPI(Pcc))
        HEM.append(nitrous_orifice.m_dot_HEM(Pcc))
        DYER.append(nitrous_orifice.m_dot_dyer(Pcc))

    plt.plot(P_cc, SPI, label='IPA')
    plt.plot(P_cc, SPIn, label='SPI nitrous')
    plt.plot(P_cc, HEM, label='HEM nitrous')
    plt.plot(P_cc, DYER, label='DYER nitrous')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
