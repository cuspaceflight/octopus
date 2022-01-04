import numpy as np
from matplotlib import pyplot as plt

from octopus import utils


def main():
    Re = np.linspace(1, 8000, 1000)

    plt.figure()

    plt.plot(Re, fd, label='fd')
    plt.plot(Re, cf, label='cf')

    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
