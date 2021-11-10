import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid


def main():
    p0 = 18e5
    p1 = 12.5e5
    T0 = 253

    nitrous = Fluid('N2O')

    nitrous.set_state(P=p0, T=T0)
    s0 = nitrous.state.smass()
    rho0 = nitrous.state.rhomass()
    x0 = nitrous.state.Q()

    nitrous.set_state(P=p1, S=s0)
    rho1 = nitrous.state.rhomass()
    x1 = nitrous.state.Q()

    D = np.linspace(rho0, rho1, 1000)[1:]
    Q = np.zeros_like(D)
    for i, d in enumerate(D):
        nitrous.set_state(D=d, S=s0)
        Q[i] = nitrous.state.Q()

    print(f'{rho0:.1f},{rho1:.1f}')
    print(f'{max(x0, 0):.3f},{x1:.3f}')

    fig,ax = plt.subplots()
    ax.plot(D, Q)
    ax.invert_xaxis()
    plt.show()

if __name__ == "__main__":
    main()
