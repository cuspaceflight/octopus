import time

from matplotlib import pyplot as plt
from numpy import linspace

from octopus import Fluid, Orifice


def main():
    fluid = Fluid('N2O', P=18e5)
    orifice = Orifice(fluid, 1e-2, 1e-3)

    P_cc = linspace(0, 19e5, 50)
    SPI = []
    HEM = []
    DYER = []
    t0 = time.time()
    for Pcc in P_cc:
        SPI.append(orifice.m_dot_SPI(Pcc))
        HEM.append(orifice.m_dot_HEM(Pcc))
        DYER.append(orifice.m_dot_dyer(Pcc))
    dt = time.time() - t0
    print(dt, len(P_cc))

    plt.plot(P_cc, SPI, label='SPI')
    plt.plot(P_cc, HEM, label='HEM')
    plt.plot(P_cc, DYER, label='DYER')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.show()


if __name__ == "__main__":
    main()
