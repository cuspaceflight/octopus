import numpy as np

from octopus import Fluid, Manifold, Nist, Orifice, PropertySource


def main():
    p0 = 18e5
    T0 = 253

    pcc = 12.5e5

    m_dot = 0.75
    OF = 3.5

    alpha = (np.pi / 180) * 20  # chosen alpha=20

    m_dot_o = m_dot * OF / (1 + OF)
    m_dot_f = m_dot * 1 / (1 + OF)

    # nitrous fluid handler
    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0, T=T0)

    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps, A=0.25 * np.pi * 25e-3 ** 2)
    ox_orifice = Orifice(manifold=ox_manifold, A=1, orifice_type=Orifice.ANNULAR)  # setup orifice with A=1 to find G

    # ipa fluid handler
    ipa_mu = 2.4e-3
    ipa_water_content = 40  # %
    data = Nist('ipa')
    concentration, density = data.get_fields('concentration', 'density')
    ipa_rho = 1000 * np.interp(100 - ipa_water_content, concentration, density)
    ipa_Cd = 0.7

    # fixed geometry definition
    inner_diam = 20e-3  # mm

    # find nitrous mass flow per unit area (assume linear scaling for inlet, little viscous effects in annulus)
    ox_G = ox_orifice.m_dot_dyer(P_cc=12.5e5)
    ox_A = m_dot_o/ox_G
    # find annular gap dimensions from flow area
    outer_diam = np.sqrt(4 * ox_A / np.pi + inner_diam ** 2)
    # find flow velocity
    ox_v = ox_G/nitrous.state.rhomass()

    # find ipa flow area directly from SPI model
    ipa_G = ipa_Cd * np.sqrt(2 * ipa_rho * (p0 - pcc))
    ipa_A = m_dot_f / ipa_G
    # find pintle total opening distance
    pintle_tod = ipa_A/(np.pi*inner_diam)
    # find flow velocity
    ipa_v = ipa_G/ipa_rho

    # find spray half angle
    theta = np.arctan(ipa_v*m_dot_f*np.cos(alpha)/(ox_v*m_dot_o+ipa_v*m_dot_f*np.sin(alpha)))

    print(f'inner diameter: {inner_diam*1000:.1f} mm\n'
          f'outer diameter: {outer_diam*1000:.1f} mm\n'
          f'annular gap: {0.5*(outer_diam-inner_diam)*1000:.2f} mm\n'
          f'total opening distance: {pintle_tod*1000:.1f} mm\n\n'
          
          f'oxidiser mass flow rate: {m_dot_o:.3f} kg/s\n'
          f'fuel mass flow rate: {m_dot_f:.3f} kg/s\n'
          f'oxidiser injection velocity: {ox_v:.1f} m/s\n'
          f'fuel injection velocity: {ipa_v:.1f} m/s\n\n'
          
          f'spray cone half angle: {180*theta/np.pi:.1f}°')


if __name__ == "__main__":
    main()
