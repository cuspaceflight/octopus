import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid, Orifice, Manifold, PropertySource


def main():
    nitrous = Fluid('N2O')
    ps = PropertySource(p=15e5, T=250)
    ox_manifold = Manifold(fluid=nitrous.name, parent=ps)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, D=1e-3, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    T = 298  # K
    psat = nitrous.psat([T])[0]
    ps._T = T

    p0 = np.linspace(psat + 100, psat + 3e5, 20)[::-1]
    P_cc = np.linspace(psat - 5e5, psat - 100, 1000)

    M = []
    for s in p0:
        ps._p = s
        Mrow = []
        for p in P_cc:
            Mrow.append(ox_orifice.m_dot(p))
        M.append(Mrow)

    plt.figure(1)
    plt.title('Dyer mass flow, increasing supercharge ^')
    [plt.plot(P_cc / 101325, m, color=((p - psat) / (p0[0] - psat), 0, 1 - (p - psat) / (p0[0] - psat)),
              label=f'supercharge: {(p - psat) / 101325:.1f} bar') for p, m in zip(p0, M)]
    plt.xlabel('Chamber Pressure/bar')
    plt.ylabel('Mass Flow/kgs$^{-1}$')
    plt.legend()

    plt.show()


if __name__ == "__main__":
    main()
