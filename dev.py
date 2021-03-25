import thermo

import octopus as octo


def main():
    # p_COM = 15e6
    # nitrous_saturated = thermo.chemical.Chemical('nitrous oxide', P=1902e3, T=255)
    # orifice = octo.Orifice(nitrous_saturated, octo.STRAIGHT, 1, 1)
    # print(orifice.C_fo,orifice.T_o,orifice.p_o,orifice.v_fo,orifice.v_fgo,orifice.h_fgo)
    # print(f'omega: {orifice.omega}, v: {orifice.v(0.85*1902e3)}')

    # h1 = octo.saturation_h(octo.NITROUS_OXIDE, 270, 'l')
    # print(f'specific enthalpy at inlet: {h1}, should be: -311')

    # For testing WIP code

    fuel_name = "isopropanol"
    oxidiser_name = "nitrous oxide"

    fuel = octo.Fluid(fuel_name, source="thermo")
    oxidiser = octo.Fluid(oxidiser_name)
    # Propellants

    water = octo.Fluid("water")
    # Additives

    fuel_mix = octo.Fluid(name_1=fuel_name, name_2="water", mf_1=0.9, mf_2=0.1)
    # Test mixture

    fuel_T = 330
    fuel_p = 18E5
    ox_T = 280
    ox_p = 40E5
    # Manifold inlet conditions, Kelvin and Pascal

    fuel_man = octo.Manifold(fluid=fuel, T_inlet=fuel_T, p_inlet=fuel_p)
    fuel_mix_man = octo.Manifold(fluid=fuel_mix, T_inlet=fuel_T, p_inlet=fuel_p)
    oxidiser_man = octo.Manifold(fluid=oxidiser, T_inlet=ox_T, p_inlet=ox_p)
    # Create (empty) manifolds for elements
    # fuel_mix_man is just to test the class properly handles a FluidMixture

    print(f'densities: IPA:{fuel_man.rho}, 90% IPA:{fuel_mix_man.rho}, nitrous oxide:{oxidiser_man.rho}')
    print(f'heat capacities: IPA:{fuel_man.Cp}, 90% IPA:{fuel_mix_man.Cp}, nitrous oxide:{oxidiser_man.Cp}')


if __name__ == "__main__":
    main()
