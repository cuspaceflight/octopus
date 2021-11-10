import numpy as np

from octopus import Fluid, Manifold, Nist, Orifice, PropertySource


def main():
    p0 = 18e5
    T0 = 253

    pcc = 12.5e5

    m_dot = 0.75
    OF = 3.5

    m_dot_o_target = m_dot * OF / (1 + OF)
    m_dot_f_target = m_dot*1/(1+OF)

    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0, T=T0)
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, A=1, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    print(f'feed pressure: {p0 / 100000} bar')
    print(f'nitrous Tsat at {p0 / 100000} bar: {nitrous.tsat([p0])[0]:.1f}\n')

    data = Nist('ipa')
    concentration, density = data.get_fields('concentration', 'density')
    ipa_density = 1000 * np.interp(80, concentration, density)
    print(f'ipa_density: {ipa_density:.1f}kg/m3\n')

    # OXIDISER
    m_flux_o = ox_orifice.m_dot_dyer(pcc)
    A_o = m_dot_o_target / m_flux_o  # A is area of cylinder to inject over
    p_o = m_dot_o_target * m_flux_o / nitrous.state.rhomass()

    r_pintle = 15 * 10 ** -3  # radius of 30mm
    # r_pintle = np.linspace(5 * 10 ** -3, 40 * 10 ** -3, 1000)  # vary radius from 2-100 mm
    h = A_o / (2 * np.pi * r_pintle)

    # FUEL
    m_flux_f = np.sqrt(2 * ipa_density * (p0 - pcc))
    A_f = m_dot_f_target / m_flux_f
    p_f = m_dot_f_target * m_flux_f / ipa_density

    r_annular = np.sqrt((A_f + np.pi * r_pintle ** 2) / np.pi)

    print(f'Oxidiser injection area: {1e6 * A_o:.0f}mm^2\n'
          f'Fuel injection area: {1e6 * A_f:.0f}mm^2\n')

    print(f'pintle opening height: {1000 * h:.2f}mm')
    print(f'annular gap: {1000 * (r_annular - r_pintle):.2f}mm\n')

    print(f'momentum ratio: {p_o:.2f}/{p_f:.2f} = {p_o / p_f:.2f}\n')

    theta = (np.pi / 180) * 20  # np.linspace(0, np.pi / 2, 100)  # chosen theta=20
    alpha = np.arctan(1 / (np.tan(theta) + p_f / (p_o * np.cos(theta))))
    print(f'pintle tip angle: {(180 / np.pi) * theta:.0f} deg\n'
          f'spray cone half angle: {(180 / np.pi) * alpha:.0f} deg')


if __name__ == "__main__":
    main()
