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
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps, A=0.25*np.pi * 25e-3 ** 2)
    gaps = np.linspace(3.2e-3, 3.5e-3, 1)
    ps = []
    plt.figure(0)

    for gap in gaps:
        r_inner = 25e-3
        annular_gap = gap  # mm
        r_outer = r_inner + annular_gap
        Dh = 2 * annular_gap
        A = np.pi * (r_outer ** 2 - r_inner ** 2)
        ox_orifice = Orifice(manifold=ox_manifold, L=25e-3, D=Dh, A=A, orifice_type=Orifice.STRAIGHT, Cd=0.7)
        orig, new = ox_orifice.p_patel(m_dot_o)
        plt.plot([cell.pos for cell in orig],[cell.v for cell in orig])
        plt.plot([cell.pos for cell in new],[cell.v for cell in new])

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
