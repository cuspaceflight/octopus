import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid, Manifold, Orifice, PropertySource


def main():
    p0 = 18e5
    T0 = 253

    pcc = 12.5e5

    m_dot = 0.75
    OF = 3.5

    alpha = (np.pi / 180) * 20  # chosen alpha=20

    m_dot_o = m_dot * OF / (1 + OF)
    m_dot_f = m_dot * 1 / (1 + OF)

    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0, T=T0)
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps, A=0.25 * np.pi * 25e-3 ** 2)
    gaps = np.linspace(3.2e-3, 3.5e-3, 1)
    ps = []
    fig, axs = plt.subplots(2,2)

    for gap in gaps:
        r_inner = 25e-3
        annular_gap = gap  # mm
        r_outer = r_inner + annular_gap
        Dh = 2 * annular_gap
        A = np.pi * (r_outer ** 2 - r_inner ** 2)
        ox_orifice = Orifice(manifold=ox_manifold, L=25e-3, D=Dh, A=A, orifice_type=Orifice.STRAIGHT, Cd=1)
        p, v, s, d = ox_orifice.p_patel_dyer(m_dot_o)

        [axs[0,0].plot(x, label=f'pressure {i}') for i, x in enumerate(p)]
        axs[0,0].set_ylabel('pressure/Pa')
        axs[0,0].legend()

        [axs[0,1].plot(x, label=f'velocity {i}') for i, x in enumerate(v)]
        axs[0,1].set_ylabel('velocity/ms-1')
        axs[0,1].legend()

        [axs[1,0].plot(x, label=f'entropy {i}') for i, x in enumerate(s)]
        axs[1,0].set_ylabel('entropy/kJ/kg.K')
        axs[1,0].legend()

        [axs[1, 1].plot(x, label=f'density {i}') for i, x in enumerate(d)]
        axs[1, 1].set_ylabel('density/kgm-3')
        axs[1, 1].legend()

    plt.xlabel('distance/m')
    plt.show()


if __name__ == "__main__":
    main()
