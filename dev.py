from matplotlib import pyplot as plt
from numpy import linspace

import octopus


def main():
    fluid = octopus.Fluid('nitrous oxide', P=50e5, T=256, method="helmholz")
    orifice_waxman = octopus.Orifice(fluid, 1e-2, 1e-3, orifice_type=1)
    orifice_straight = octopus.Orifice(fluid, 1e-2, 1e-3, orifice_type=0)

    # Used to check against the injector Waxman built
    fluid2 = octopus.Fluid('nitrous oxide', P=48.4e5, T=256, method="helmholz")
    check_orifice = octopus.Orifice(fluid2, 1e-2, 2.15e-3, orifice_type=1)

    [print(f'{key}: {val}') for key, val in check_orifice.m_dot_waxman(37.1e5).items()]
    print(f'D2: {check_orifice.D}')

    P_cc = linspace(0, 45e5, 50)
    SPI = []
    HEM = []
    DYER = []
    for Pcc in P_cc:
        SPI.append(orifice_straight.m_dot_SPI(Pcc))
        HEM.append(orifice_straight.m_dot_HEM(Pcc))
        DYER.append(orifice_straight.m_dot_dyer(Pcc))

    plt.plot(P_cc, SPI, label='SPI')
    plt.plot(P_cc, HEM, label='HEM')
    plt.plot(P_cc, DYER, label='DYER')
    plt.xlabel('Downstream pressure (Pa)')
    plt.ylabel('Mass Flow Rate (kg/s)')
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
