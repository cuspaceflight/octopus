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
    fluid = Fluid('nitrous oxide', P=18e5)
    orifice = Orifice(fluid, 1e-2, 1e-3)

    P_cc = linspace(0, 19e5, 50)
    SPI = []
    HEM = []
    DYER = []
    for Pcc in P_cc:
        SPI.append(orifice.m_dot_SPI(Pcc))
        HEM.append(orifice.m_dot_HEM(Pcc))
        DYER.append(orifice.m_dot_dyer(Pcc))

    plt.plot(P_cc, SPI, label='SPI')
    plt.plot(P_cc, HEM, label='HEM')
    plt.plot(P_cc, DYER, label='DYER')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.show()


if __name__ == "__main__":
    main()
