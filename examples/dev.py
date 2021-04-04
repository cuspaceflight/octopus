from octopus import Fluid, Orifice

"""LEGACY CPL CODE

if not 183.15 <= T <= 303.15:
    raise ValueError(f"Temperature ({T} K) out of range")
Tr = 309.57  # Find the reduced temperature (T / T_crit)
return 2.49973 * (1 + 0.023454 / (1 - Tr) - 3.80136 * (1 - Tr) +
                  13.0945 * (1 - Tr) ** 2 - 14.5180 * (1 - Tr) ** 3)
"""


def main():
    nitrous = Fluid('nitrous oxide', P=18e5, method='helmholz')
    isopropanol = Fluid('isopropanol', P=18e5, method='thermo')
    nitrous_orifice = Orifice(nitrous, 1e-2, 1e-3)
    ipa_orifice = Orifice(isopropanol, 1e-2, 1e-3)

    T = 193
    rho = nitrous.rho_g(T)

    print(nitrous.cp(rho, T) - nitrous.cv(rho, T),nitrous.MW,nitrous.R_specific)
    print(nitrous.cp(rho, T) / nitrous.cv(rho, T))


if __name__ == "__main__":
    main()
