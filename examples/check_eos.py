import numpy as np
from matplotlib import pyplot as plt

from octopus import Fluid


def main():
    nitrous_heos = Fluid('N2O', 'HEOS')
    nitrous_pr = Fluid('N2O', 'PR')

    T = np.linspace(182.23, 300, 1000)
    ad = 0.5*(T-270)

    n_heos_rhog = nitrous_heos.rhog(T)
    n_heos_rhol = nitrous_heos.rhol(T)

    n_pr_rhog = np.array(nitrous_pr.rhog(T))
    n_pr_rhol = np.array(nitrous_pr.rhol(T))+ad

    fig, ax = plt.subplots(1)

    ax.plot(T, n_heos_rhog, color='green', label='heos rhog')
    ax.plot(T, n_pr_rhog, color='blue', label='pr rhog')

    ax.plot(T, n_heos_rhol, color='green', label='heos rhol')
    ax.plot(T, n_pr_rhol, color='blue', label='pr rhol')

    error_l = 100 * abs((n_heos_rhol - n_pr_rhol) / n_heos_rhol)
    error_g = 100 * abs((n_heos_rhog - n_pr_rhog) / n_heos_rhog)
    print('density: error = (HELMHOLZ-PR)/HELMHOLZ')
    print(f'saturated liquid error:'
          f' mean = {np.mean(error_l):.1f}%,'
          f' max = {np.max(error_l):.1f}%,'
          f' at 250K = {np.interp(250, T, error_l):.1f}%')
    print(f'saturated vapour error: '
          f' mean = {np.mean(error_g):.1f}%,'
          f' max = {np.max(error_g):.1f}%,'
          f' at 250K = {np.interp(250, T, error_g):.1f}%')

    ax.set_xlabel('temperature/K')
    ax.set_ylabel('density/kgm-3')

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
