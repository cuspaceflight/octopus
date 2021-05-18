import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid


def main():
    nitrous = Fluid('N2O')
    T = np.linspace(150, 310, 1000)

    rhol = nitrous.rhol(T)
    rhog = nitrous.rhog(T)

    hl = nitrous.hl(T)
    hg = nitrous.hg(T)

    sl = nitrous.sl(T)
    sg = nitrous.sg(T)

    psat = nitrous.psat(T)

    plt.figure(1)
    plt.title('Saturation Density')
    plt.plot(T, rhol, color='blue', label='liquid')
    plt.plot(T, rhog, color='red', label='vapour')
    plt.legend()

    plt.figure(2)
    plt.title('Saturation Enthalpy')
    plt.plot(T, hl, color='blue', label='liquid')
    plt.plot(T, hg, color='red', label='vapour')
    plt.legend()

    plt.figure(3)
    plt.title('Saturation Entropy')
    plt.plot(T, sl, color='blue', label='liquid')
    plt.plot(T, sg, color='red', label='vapour')
    plt.legend()

    plt.figure(4)
    plt.title('Saturation Pressure')
    plt.plot(T, psat, color='blue')

    plt.show()


if __name__ == "__main__":
    main()
