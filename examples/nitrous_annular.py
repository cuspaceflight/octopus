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
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps)
    r_inner = 25e-3
    annular_gap = 3.05e-3 # mm
    r_outer = r_inner + annular_gap
    Dh = 2 * annular_gap
    A = np.pi * (r_outer ** 2 - r_inner ** 2)
    print(A)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, D=Dh, A=A, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    plt.figure(0)
    plt.plot(ox_orifice.dP_HEM_frictional(m_dot_o))
    plt.show()


if __name__ == "__main__":
    main()
