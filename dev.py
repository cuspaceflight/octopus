from matplotlib import pyplot as plt
from numpy import array
from numpy import linspace

from octopus import Fluid

fluid = Fluid('N2O', T=250, P=18e5)
fluid.calculate()


def main():
    T = array(linspace(260, 309, 100))
    r_helm = [fluid.rho_sat(t) for t in T]
    rg_helm = [rg for rg, rl in r_helm]
    rl_helm = [rl for rg, rl in r_helm]
    rg_thermo = [fluid.rho_g(t) for t in T]
    rl_thermo = [fluid.rho_l(t) for t in T]

    plt.plot(T, rg_helm, color='red')
    plt.plot(T, rl_helm, color='red')
    plt.plot(T, rg_thermo, color='green')
    plt.plot(T, rl_thermo, color='green')
    plt.scatter(fluid.Tc, fluid.rhoc, color='blue')
    plt.show()


if __name__ == "__main__":
    main()
