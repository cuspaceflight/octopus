import octopus as octo

# For testing WIP code

fuel_name = "isopropanol"
oxidiser_name = "nitrous oxide"

fuel = octo.Fluid(fuel_name)
oxidiser = octo.Fluid(fuel_name)
# Propellants

water = octo.Fluid("water")
# Additives

fuel_mix = octo.FluidMixture("water", fuel_name, 0.1, 0.9,)

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
