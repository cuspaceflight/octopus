import matplotlib.pyplot as plt
import numpy as np

import octopus as octo


def main():
    # Configuration of the fluid and orifice

    T_N2O = 257  # K

    nitrous = octo.Fluid('nitrous oxide', P=25e5, T=T_N2O, method='helmholz')

    test_orifice = octo.Orifice(nitrous, 1e-2, 1e-3, orifice_type=0)

    # Nitrous oxide is injected into the lower chamber, initially at 1 bar, no nitrous oxide.
    # The introduction of nitrous is assumed to be quasi-equlibrium and adiabatic (so isentropic).
    # Both air and nitrous are assumed to behave ideally in the lower chamber.

    # State 0: p0 = 1 bar, m_N2O = 0.

    R = 8.31446  # J/mol/K
    M_N2O = 0.04401  # Kg/mol

    D_in = 0.05  # m
    L = 1  # m
    P0 = 1E5  # Pa
    P = [P0]
    V_chamber = np.pi * L * D_in ** 2 / 4

    R_N2O = R / M_N2O

    # For each timestep, the downstream pressure is reevlauated for the next mass flow rate calculation
    # The nitrous is assumed to quickly reach a pure vapour form in the lower chamber, with the same
    # temperature at which it entered the injector (T_N2O).

    P_max = 20E5
    m_N2O = [0]
    t = [0]
    dt = 0.01

    i = 0
    j = 0

    while P[-1] < P_max:
        # Update time and iteration record
        i += 1
        t.append(t[-1] + dt)

        # Calculate the mass flow increment using the most recent backpressure
        # and add it to the total nitrous mass in the lower chamber

        # Check if Dyer is having problems
        if test_orifice.m_dot_dyer(P[-1]) is not None:
            dm = dt * test_orifice.m_dot_dyer(P[-1])
        else:
            dm = m_N2O[-1] - m_N2O[-2]
            print(f"Dyer did not return a valid result for iteration {i}")
            j += 1
            # If so, reuse the last dm and hope it fixes itself next time
        m_N2O.append(m_N2O[-1] + dm)

        # Dalton's law of partial pressures
        P.append(P0 + m_N2O[-1] * R_N2O * T_N2O / V_chamber)

        if i % 100 == 0:
            print(f"Iteration {i}, t = {t[-1]:.2f} s, last nitrous mass increase {dm:.7f} kg")
            print(f"Total nitrous mass {m_N2O[-1]:.4f} kg, current backpressure {P[-1] / 1E5:.4f} bar")
            print("")

    print(f"Iteration {i}, last nitrous mass increase {dm:.7f} kg")
    print(f"Total nitrous mass {m_N2O[-1]:.4f} kg, current backpressure {P[-1] / 1E5:.4f} bar")
    print(f"Backpressure sweep time: {t[-1]:.4f} s, Dyer failed {j} of {len(t)} times")

    plt.plot(t, [p / 1E5 for p in P])
    plt.xlabel("Time, s")
    plt.ylabel("Backpressure, bar")
    plt.title("Nitrous cold flow lower chamber backpressure sweep")
    plt.show()


if __name__ == "__main__":
    main()
