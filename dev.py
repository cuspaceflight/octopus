import octopus as octo

# For testing WIP code

fuel_name = "isopropanol"
oxidiser_name = "nitrous oxide"

fuel = octo.Fluid(fuel_name)
oxidiser = octo.Fluid(fuel_name)

fuel_mix = octo.FluidMixture(fuel_name, "water", 0.9, 0.1)
