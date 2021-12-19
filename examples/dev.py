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

    pcc = 12.5e5

    m_dot = 0.75
    OF = 3.5

    alpha = (np.pi / 180) * 20  # chosen alpha=20

    m_dot_o_target = m_dot * OF / (1 + OF)
    m_dot_f_target = m_dot * 1 / (1 + OF)

    print(f'mass flows:\n\toxidiser - {m_dot_o_target:.3f}\n\tfuel - {m_dot_f_target:.3f}')

    nitrous = Fluid('N2O')
    nitrous_ps = PropertySource(p=p0, T=T0)
    nitrous.set_state(P=p0,T=T0)
    ox_manifold = Manifold(fluid=nitrous, parent=nitrous_ps)
    ox_orifice = Orifice(manifold=ox_manifold, L=1e-2, A=1, orifice_type=Orifice.STRAIGHT, Cd=0.7)

    print(f'feed pressure: {p0 / 100000} bar')
    print(f'nitrous Tsat at {p0 / 100000} bar: {nitrous.tsat([p0])[0]:.1f}\n')

    data = Nist('ipa')
    concentration, density = data.get_fields('concentration', 'density')
    ipa_density = 1000 * np.interp(80, concentration, density)
    ipa_mu = 2.4e-3
    print(f'ipa_density: {ipa_density:.1f}kg/m3\n')

    # OXIDISER
    print(f'mu before: {nitrous.viscosity()}')
    m_flux_o = ox_orifice.m_dot_dyer(pcc)
    A_o = m_dot_o_target / m_flux_o  # A is area of cylinder to inject over
    p_o = m_dot_o_target * m_flux_o / nitrous.state.rhomass()
    V_o = m_dot_o_target / (nitrous.state.rhomass() * A_o)
    print(f'mu after: {nitrous.viscosity()}')

    r_pintle = 10e-3  # radius of 20mm
    L = 20e-3  # length of annular gap
    h = A_o / (2 * np.pi * r_pintle * np.cos(alpha))

    # FUEL
    res = scipy.optimize.root_scalar(f=dp_annular_gap,
                                     args=(r_pintle, m_dot_f_target, L, ipa_density, ipa_mu, p0 - pcc),
                                     x0=r_pintle * 1.01, x1=r_pintle * 1.02)
    print(res)
    r_annular = res.root

    A_f = np.pi * (r_annular ** 2 - r_pintle ** 2)  # flow area
    D = 2 * (r_annular - r_pintle)
    m_flux_f = m_dot_f_target / A_f
    p_f = m_dot_f_target * m_flux_f / ipa_density
    V_f = m_dot_f_target / (ipa_density * A_f)

    Re_f = ipa_density * V_f * D / ipa_mu  # Reynolds number

    print(f'Oxidiser injection area: {1e6 * A_o:.2f}mm^2\n'
          f'Fuel injection area: {1e6 * A_f:.2f}mm^2\n')

    print(f'pintle opening height: {1000 * h:.4f}mm')
    print(f'annular gap: {1000 * (r_annular - r_pintle):.4f}mm\n')

    print(f'Oxidiser injection velocity: {V_o:.1f}m/s')
    print(f'Fuel injection velocity: {V_f:.1f}m/s\n')

    print(f'Fuel Reynolds number at exit: {Re_f:.0f}\n')

    print(f'momentum ratio: {p_o:.2f}/{p_f:.2f} = {p_o / p_f:.2f}\n')

    theta = np.arctan(1 / (np.tan(alpha) + p_f / (p_o * np.cos(alpha))))
    print(f'pintle tip angle: {(180 / np.pi) * alpha:.0f} deg\n'
          f'spray cone half angle: {(180 / np.pi) * theta:.0f} deg')


if __name__ == "__main__":
    main()
