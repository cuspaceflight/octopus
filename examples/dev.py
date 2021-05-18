import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid, Orifice, Manifold, PropertySource


def main():
    nitrous = Fluid('N2O')
    ps = PropertySource(p=15e5, T=250)
    ox_manifold = Manifold(fluid=nitrous.name, parent=ps)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, D=1e-3, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    P1 = np.linspace(1e5, 15e5, 1000)
    P0 = np.linspace(15e5, 20e5, 10)
    M = []

    for p0 in P0:
        ps._p = p0
        Mrow = []
        for p1 in P1:
            Mrow.append(ox_orifice.m_dot(P_cc=p1))
        M.append(Mrow)

    T = np.linspace(150, 310, 1000)

    rhol = nitrous.rhol(T)
    rhog = nitrous.rhog(T)
    psat = nitrous.psat(T)

    plt.figure(1)
    plt.title('Dyer mass flow, increasing supercharge ^')
    [plt.plot(P1, m, color='black') for m, p0 in zip(M, P0)] #label=f'p0: {p0}'
    plt.legend()
    '''
    plt.figure(2)
    plt.title('Density')
    plt.plot(T, rhol, color='blue')
    plt.plot(T, rhog, color='red')

    plt.figure(3)
    plt.title('Pressure')
    plt.plot(T, psat, color='blue')
    '''
    plt.show()


if __name__ == "__main__":
    main()
