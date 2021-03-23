import matplotlib.pyplot as plt
import numpy as np


def mdot_spi(Cd, A, n, rho, p_chamber):
    # Calculates and plots a single phase incompressible (SPI) approximation of the mass flow.

    N = 1000
    dp = np.linspace(0.15 * p_chamber, 0.2 * p_chamber, N)  # 25-20% of chamber pressure

    mdot = n * Cd * A * np.sqrt(2 * rho * dp)  # aiming for 1.3 kg/s for fuel

    ax = plt.subplot(1, 1, 1)
    ax.plot(dp, mdot)
    ax.set_xlabel("Pressure Drop [Pa]")
    ax.set_ylabel("Mass Flow Rate [kgs$^{-1}$]")
    ax.set_title("Fuel Flow Rate (SPI)")

    plt.show()


def mdot_hem(A):
    # produces mass flow given a total area, assuming choked flow in HEM

    h_fgo = 1000 * ((-69.8) - (-355))  # J/kg
    v_fgo = (1 / 40.11) - (1 / 1014.8)  # m3/kg
    C_fo = 1915  # J/Kg.K
    T_o = 273 - 25  # K
    return (h_fgo / v_fgo) * np.power(C_fo * T_o, -0.5) * A


def req_A(mdot):
    # produces the required total area for a given mass flow using HEM

    h_fgo = 1000 * ((-69.3) - (-345))  # J/kg
    v_fgo = (1 / 46.82) - (1 / 995.4)  # m3/kg
    C_fo = 1957  # J/Kg.K
    T_o = 273 - 20  # K
    return (v_fgo / h_fgo) * np.power(C_fo * T_o, 0.5) * mdot


def delta_p(mdot, A, n, Cd, rho):
    return np.power(mdot / (A * n * Cd), 2) / (2 * rho)


if __name__ == "__main__":
    Cd = 0.7  # sensible discharge coefficient
    A = np.pi * 0.5e-3 ** 2  # area of a 0.5mm radius orifice
    rho = 995.4  # kg/m3
    mdot=4.1
    p_chamber = 1500000  # 1.5MPa, 15bar

    # mdot_spi(Cd, A, n, rho, p_chamber)
    calc_A = req_A(mdot)
    n=int(calc_A/A)
    dp=delta_p(mdot,A,n,Cd,rho)
    print(f'required area for {mdot}kg/s: {round(calc_A,7)}m^2. Assuming 1mm diameter, {n} orifices required')
    print(f'required pressure drop for {mdot}kg/s over {n} {round(A,7)}m^2 orifices: {round(dp,4)}Pa')
