import numpy as np
from matplotlib import pyplot as plt

from octopus import utils


def main():
    Re = np.linspace(500, 8000, 1000)

    plt.figure()

    plt.plot(Re, [utils.fd_solved(r) for r in Re], label='fd_solved')
    plt.plot(Re, [utils.fd(r) for r in Re], label='fd')

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
