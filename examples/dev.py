import numpy as np
import scipy.optimize
import scipy.special

from octopus import Fluid, Manifold, Nist, Orifice, PropertySource


def dp_annular_gap(r_outer, r_inner, mdot, L, rho, mu, dp=0):
    A = np.pi * (r_outer ** 2 - r_inner ** 2)
    D = 2 * (r_outer - r_inner)
    V = mdot / (rho * A)

    Re = rho * V * D / mu
    return 0.5 * rho * (1 + fd(Re) * L / D) * V ** 2 - dp


# often fd/4 is used in literature: check
def fd(Re):
    laminar = 64 / Re
    turbulent_smooth = 1 / np.real(0.838 * scipy.special.lambertw(0.629 * Re)) ** 2
    if Re < 2000:
        return laminar
    elif Re > 4000:
        return turbulent_smooth
    else:
        return (laminar * (4000 - Re) + turbulent_smooth * (Re - 2000)) / (4000 - 2000)


def main():
    roughness = 5e-6
    p0 = 18e5
    T0 = 253

    pcc = 10e5

    m_dot = 0.803
    OF = 3.5

    alpha = (np.pi / 180) * 20  # chosen alpha=20

    m_dot_o = m_dot * OF / (1 + OF)
    m_dot_f = m_dot * 1 / (1 + OF)

    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0, T=T0)
    nitrous.set_state(P=p0, T=T0)
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps, A=1)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, A=1, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    data = Nist('ipa')
    concentration, density = data.get_fields('concentration', 'density')
    ipa_rho = 1000 * np.interp(80, concentration, density)
    ipa_mu = 2.4e-3

    # OXIDISER
    print(nitrous.state.rhomass())
    print(nitrous.psat([253])[0])
    ox_G = ox_orifice.m_dot_dyer(pcc)
    ox_A = m_dot_o / ox_G  # A is area of cylinder to inject over
    ox_v = m_dot_o / (nitrous.state.rhomass() * ox_A)


    r_pintle = 10e-3  # diameter of 20mm
    L = 10e-3  # length of annular gap
    h = ox_A / (2 * np.pi * r_pintle * np.cos(alpha))

    # FUEL
    res = scipy.optimize.root_scalar(f=dp_annular_gap,
                                     args=(r_pintle, m_dot_f, L, ipa_rho, ipa_mu, p0 - pcc),
                                     x0=r_pintle * 1.01,
                                     x1=r_pintle * 1.02)

    r_annular = res.root

    ipa_A = np.pi * (r_annular ** 2 - r_pintle ** 2)  # flow area
    D = 2 * (r_annular - r_pintle)
    ipa_v = m_dot_f / (ipa_rho * ipa_A)

    TMR = ox_v * m_dot_o * np.cos(alpha) / (ipa_v * m_dot_f + ox_v * m_dot_o * np.sin(alpha))
    theta = np.arctan(TMR)

    print(f'inner diameter: {2 * r_pintle * 1000:.2f} mm\n'
          f'outer diameter: {2 * r_annular * 1000:.2f} mm\n'
          f'annular gap: {(r_annular - r_pintle) * 1000:.2f} mm\n'
          f'total opening distance: {h * 1000:.2f} mm\n\n'

          f'oxidiser injection area: {1e6*ox_A:.2f}mm^2\n'
          f'fuel injection area: {1e6*ipa_A:.2f}mm^2\n\n'

          f'oxidiser mass flow rate: {m_dot_o:.3f} kg/s\n'
          f'fuel mass flow rate: {m_dot_f:.3f} kg/s\n'
          f'oxidiser injection velocity: {ox_v:.1f} m/s\n'
          f'fuel injection velocity: {ipa_v:.1f} m/s\n'
          f'TMR: {TMR:.2f}\n\n'

          f'spray cone half angle: {180 * theta / np.pi:.1f}°')


if __name__ == "__main__":
    main()
