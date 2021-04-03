from matplotlib import pyplot as plt
from numpy import linspace

from octopus import Fluid, Orifice, STRAIGHT


def main():
    fluid = Fluid('N2O', T=250, P=18e5)
    orifice = Orifice(fluid, 0, 15e5, STRAIGHT, 1e-2, 1e-3)

    P_cc = linspace(0, 19e5, 100)
    SPI = []
    HEM = []
    DYER = []
    for Pcc in P_cc:
        orifice.P_cc = Pcc
        SPI.append(orifice.m_dot_SPI())
        HEM.append(orifice.m_dot_HEM())
        DYER.append(orifice.m_dot_dyer())

    plt.plot(P_cc, SPI, label='SPI')
    plt.plot(P_cc, HEM, label='HEM')
    plt.plot(P_cc, DYER, label='DYER')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.show()


if __name__ == "__main__":
    main()
