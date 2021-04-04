from matplotlib import pyplot as plt
from numpy import linspace

from octopus import Fluid, Orifice

"""LEGACY CPL CODE

if not 183.15 <= T <= 303.15:
    raise ValueError(f"Temperature ({T} K) out of range")
Tr = 309.57  # Find the reduced temperature (T / T_crit)
return 2.49973 * (1 + 0.023454 / (1 - Tr) - 3.80136 * (1 - Tr) +
                  13.0945 * (1 - Tr) ** 2 - 14.5180 * (1 - Tr) ** 3)
"""


def main():
    nitrous = Fluid('nitrous oxide', P=18e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, method='thermo')
    nitrous_orifice = Orifice(nitrous, 1e-2, 1e-3)
    ipa_orifice = Orifice(isopropanol, 1e-2, 1e-3)

    P_cc = linspace(0, 19e5, 50)
    SPI = []
    DYER = []
    for Pcc in P_cc:
        SPI.append(ipa_orifice.m_dot_SPI(Pcc))
        DYER.append(nitrous_orifice.m_dot_dyer(Pcc))

    plt.plot(P_cc, SPI, label='IPA')
    plt.plot(P_cc, DYER, label='N2O')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
