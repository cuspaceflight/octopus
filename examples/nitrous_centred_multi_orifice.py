import numpy as np
import scipy.optimize
import scipy.special

from octopus import Fluid, Manifold, Nist, Orifice, PropertySource, utils


def main():
    # geometry and boundary condition set
    p0_nitrous = 27.8e5  # 17.9bar vapour pressure
    p0_ipa = p0_nitrous  # set *total* pressure of inlet fluids
    T0 = 273-20  # inlet temperature of N2O
    pcc = 10e5  # chamber pressure
    m_dot = 0.803  # total mass flow rate
    OF = 3.5  # Oxidiser/Fuel ratio
    water_content = 25  # percentage water content of fuel (Isp max at 25%)
    alpha = (np.pi / 180) * 40

    # oxidiser and fuel mass flow rates
    m_dot_o = m_dot * OF / (1 + OF)
    m_dot_f = m_dot * 1 / (1 + OF)

    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0_nitrous, T=T0)
    nitrous.set_state(P=p0_nitrous, T=T0)
    ox_rho0 = nitrous.state.rhomass()
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps, A=1)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, A=1, orifice_type=Orifice.STRAIGHT, Cd=0.47)

    data = Nist('ipa')
    concentration, density = data.get_fields('concentration', 'density')
    ipa_rho = 1000 * np.interp(100 - water_content, concentration, density)
    ipa_mu = 2.4e-3

    # OXIDISER
    ox_G = ox_orifice.m_dot_dyer(pcc)
    print(f'real ox mdot: {ox_G*31.45e-6}')
    ox_A = m_dot_o / ox_G  # A is area of cylinder to inject over
    ox_v = m_dot_o / (nitrous.state.rhomass() * ox_A)

    v_crit = np.sqrt(2 * (p0_nitrous - nitrous.psat([T0])[0]) / ox_rho0)
    A_crit = m_dot_o / (ox_rho0 * v_crit)
    D_crit = np.sqrt(4 * A_crit / np.pi)

    D_pintle = 11e-3
    sleeve_thickness = 1e-3
    slot_width = 1e-3
    num_slots = np.floor(np.pi * (D_pintle - 2 * sleeve_thickness) / slot_width / 2)
    slot_height = ox_A / num_slots / slot_width + sleeve_thickness * np.sin(alpha)

    # FUEL
    L = 30e-3
    # res = scipy.optimize.root_scalar(f=utils.dp_annular_gap,
    #                                  args=(D_pintle, m_dot_f, L, ipa_rho, ipa_mu, p0_ipa - pcc),
    #                                  x0=D_pintle * 1.01,
    #                                  x1=D_pintle * 1.02)

    D_outer = 11.94e-3  # res.root

    res = scipy.optimize.root_scalar(f=utils.dp_ipa_throttle,
                                     args=(D_outer, D_pintle, m_dot_f, ipa_rho, p0_nitrous - pcc, 0.65),
                                     x0=0.4,
                                     x1=0.41)

    rat_A = res.root

    ipa_A = np.pi * (D_outer ** 2 - D_pintle ** 2) / 4  # flow area
    D = D_outer - D_pintle
    ipa_v = m_dot_f / (ipa_rho * ipa_A)

    TMR = ox_v * m_dot_o * np.cos(alpha) / (ipa_v * m_dot_f + ox_v * m_dot_o * np.sin(alpha))
    theta = np.arctan(TMR)

    # dP_ipa_manifold = p0_nitrous - p0_ipa
    # rat_A = 1 / np.sqrt(ipa_A ** 2 * 2 * ipa_rho * dP_ipa_manifold / (m_dot_f ** 2) + 1)
    mill_sizes = 0.001 * np.arange(1.5, 4.5, 0.5)
    num_holes = np.pi * rat_A * D_pintle / mill_sizes
    for m, n in zip(mill_sizes, num_holes): print(m, n)

    print(f'OXIDISER CENTRED, MULTI-ORIFICE PINTLE\n'
          f'======================================')

    print(f'nitrous manifold pressure: {p0_nitrous / 100000:.1f}bar\n'
          f'nitrous manifold temp: {T0:.0f}K\n'
          f'ipa manifold pressure: {p0_ipa / 100000:.1f}bar\n'
          f'chamber pressure: {pcc / 100000:.1f}bar\n')

    print(f'inner diameter: {D_pintle * 1000:.2f} mm\n'
          f'outer diameter: {D_outer * 1000:.2f} mm\n'
          f'annular gap: {0.5 * (D_outer - D_pintle) * 1000:.2f} mm\n')

    print(f'Nitrous slots:\n'
          f'  dimensions: {slot_width * 1000:.1f}x{slot_height * 1000:.2f}mm\n'
          f'  count: {int(num_slots)}')

    print(f'\nFuel area reduction ratio: {rat_A:.2f}\n')
    print(f'oxidiser injection area: {1e6 * ox_A:.2f}mm^2\n'
          f'fuel injection area: {1e6 * ipa_A:.2f}mm^2\n'
          f'oxidiser critical area: {A_crit * 1e6:.2f}, diameter: {D_crit * 1000:.2f}mm\n\n'

          f'oxidiser mass flow rate: {m_dot_o:.3f} kg/s\n'
          f'fuel mass flow rate: {m_dot_f:.3f} kg/s\n'
          f'-> water content: {water_content}%\n'
          f'oxidiser injection velocity: {ox_v:.1f} m/s\n'
          f'fuel injection velocity: {ipa_v:.1f} m/s\n'
          f'TMR: {TMR:.2f}\n\n'

          f'spray cone half angle: {180 * theta / np.pi:.1f}°')


if __name__ == "__main__":
    main()
