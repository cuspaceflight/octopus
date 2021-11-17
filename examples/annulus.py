import numpy as np
import scipy.special


def mdot_annular_gap(mdot, dp, r_inner, r_outer, L, density, mu):
    A = np.pi * (r_outer ** 2 - r_inner ** 2)
    D = 2 * (r_outer - r_inner)
    V = mdot / (density * A)

    Re = density * V * D / mu
    m_flux = np.sqrt((2 * density * dp) / (1 + fd(Re) * L / D))

    return A * m_flux


def fd(Re):
    if Re < 2000:
        return 64 / Re
    else:
        return 1 / np.real(0.838 * scipy.special.lambertw(0.629 * Re)) ** 2
        # https://en.wikipedia.org/wiki/Darcy%E2%80%93Weisbach_equation?section=9#Smooth-pipe_regime


def main():
    p0 = 18e5
    pcc = 12.5e5
    m_dot = 0.75
    OF = 3.5
    m_dot_f_target = m_dot * 1 / (1 + OF)

    ipa_density = 827.3 # calculated from charts
    ipa_mu = 2.4e-3 # set data

    r_pintle = 10e-3  # radius of 20mm
    r_annular = r_pintle + 1e-3  # initial guess of 1mm gap
    L = 20e-3  # length of annular gap 20mm

    m = m_dot_f_target
    i = 0
    while True:
        i += 1
        m = mdot_annular_gap(m, p0 - pcc, r_pintle, r_annular, L, ipa_density, ipa_mu)
        r_annular += 0.0001 * (m_dot_f_target - m)
        print(m)
        if abs(m_dot_f_target - m) < 0.00001:
            break
        elif i > 100:
            raise RecursionError
    A = np.pi * (r_annular ** 2 - r_pintle ** 2)
    V = m_dot_f_target / (ipa_density * A)
    D = 2 * (r_annular - r_pintle)
    Re = ipa_density * V * D / ipa_mu  # Reynolds number

    print(f'Fuel injection area: {1e6 * A:.2f}mm^2\n')

    print(f'annular gap: {1000 * (r_annular - r_pintle):.4f}mm\n')

    print(f'Fuel injection velocity: {V:.1f}m/s')
    print(f'Fuel Reynolds number at exit: {Re:.0f}')


if __name__ == "__main__":
    main()
