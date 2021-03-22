import matplotlib.pyplot as plt
import numpy as np


def fuel_mass_flow():
    N = 1000

    Cd = 0.77                                               # sensible discharge coefficient
    A = np.pi * 0.5e-3 ** 2                          # area of a 0.5mm radius orifice
    n = 95                                            # number of orifices
    rho = 1000                                              # kg/m3
    p_chamber = 1500000                                     # 1.5MPa, 15bar
    dp = np.linspace(0.15 * p_chamber, 0.2 * p_chamber, N)  # 25-20% of chamber pressure

    mdot = n * Cd * A * np.sqrt(2 * rho * dp)                   # aiming for 1.3 kg/s for fuel

    ax = plt.subplot(1, 1, 1)
    ax.plot(dp, mdot)
    ax.set_xlabel("Pressure Drop [Pa]")
    ax.set_ylabel("Mass Flow Rate [kgs$^{-1}$]")
    ax.set_title("Fuel Flow Rate (Standard Model)")

    plt.show()


if __name__ == "__main__":
    fuel_mass_flow()
