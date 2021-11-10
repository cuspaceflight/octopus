import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid


def main():
    nitrous = Fluid('N2O')
    T = np.linspace(182.23, 309.5, 1000)
    # T = np.linspace(64, 126, 1000)

    rhol = nitrous.rhol(T)
    rhog = nitrous.rhog(T)

    hl = nitrous.hl(T)
    hg = nitrous.hg(T)

    sl = nitrous.sl(T)
    sg = nitrous.sg(T)

    psat = nitrous.psat(T)
    psat = [p * 0.00001 if p else None for p in psat]

    fig, axs = plt.subplots(2, 2, sharex=True)

    fig.suptitle('N$_2$O Saturation Data')

    axs[0, 0].plot(T, rhol, color='blue', label='liquid')
    axs[0, 0].plot(T, rhog, color='red', label='vapour')
    axs[0, 0].legend()
    axs[0, 0].set_title('Density (kg/m3)')

    axs[1, 0].plot(T, hl, color='blue', label='liquid')
    axs[1, 0].plot(T, hg, color='red', label='vapour')
    axs[1, 0].legend()
    axs[1, 0].set_title('Enthalpy (kJ/kg)')

    axs[0, 1].plot(T, sl, color='blue', label='liquid')
    axs[0, 1].plot(T, sg, color='red', label='vapour')
    axs[0, 1].legend()
    axs[0, 1].set_title('Entropy (kJ/kg.K)')

    axs[1, 1].set_title('Pressure (bar)')
    axs[1, 1].plot(T, psat, color='blue')

    plt.show()


if __name__ == "__main__":
    main()
