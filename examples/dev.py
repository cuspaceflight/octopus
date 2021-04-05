import numpy as np
from octopus import Fluid, Orifice
import matplotlib.pyplot as plt


def main():
    nitrous = Fluid('nitrous oxide', P=18e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, method='thermo')
    nitrous_orifice = Orifice(nitrous, 1e-2, 1e-3)
    ipa_orifice = Orifice(isopropanol, 1e-2, 1e-3)

    print(nitrous.Tc)

    # Used to check against the injector Waxman built
    fluid2 = Fluid('nitrous oxide', P=48.4e5, T=256, method="helmholz")
    check_orifice = Orifice(fluid2, 1e-2, 2.15e-3, orifice_type=1)

    [print(f'{key}: {val}') for key, val in check_orifice.m_dot_waxman(37.1e5).items()]
    print(f'D2: {check_orifice.D}')

    P_cc = np.linspace(0, 19e5, 50)
    SPI = []
    SPIn=[]
    HEM=[]
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
