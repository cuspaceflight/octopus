"""Implementation of data and mathematical utilities used in the `octopus.main` module"""
from collections import OrderedDict
from csv import reader
from typing import Any, List, Union

import numpy as np
import scipy

from numpy import array, nan_to_num, ndarray, real
from scipy.special import lambertw


def derivative(f: callable, axis: int, *args: Any, dx: float = 0.01) -> float:
    """Calculate the  numerical derivative of a scalar funtion f(args) with respect to args[axis].

    :param f: vector function to find derivative of with the signature f(args) -> iterable
    :param axis: index of the parameter against which the derivative of f is to be taken
    :param args: parameters to pass to f, including the one against which the derivative is to be taken
    :param dx: differentiation step

    :return: array of floats representing vector derivative f with respect to args[axis]

    """
    if axis >= len(args):
        raise ValueError("axis index too high")
    args_l = list(args)
    args_h = list(args)
    args_l[axis] -= dx
    args_h[axis] += dx

    return nan_to_num(array((f(*args_h) - f(*args_l)) / (2 * dx)))


class Nist:
    """Collection of tools for retrieving data from a NIST csv file."""

    def __init__(self, filename: str):
        """Retrieve thermodynamic data from a tab-delimited csv file.

        :param filename: filename, not including file extension, of a tab-delimited NIST datafile

        """
        with open(f'{filename}.csv', newline='') as f:
            csvreader = reader(f, delimiter='\t')
            rows = [row for row in csvreader]
            self.data = OrderedDict()
            for i, header in enumerate(rows[0]):
                self.data[header] = [float(row[i]) if row[i] != 'undefined' else None for row in rows[1:]]

    def list_fields(self) -> List[str]:
        """List fields in NIST datafile

        :returns: list of column headers in the NIST datafile

        """
        return [key for key in self.data.keys()]

    def get_fields(self, *fields: Union[str, float]) -> List[ndarray]:
        """Retrieve data under specified fields

        :param fields: titles or indices of column headers in the datafile for which data is to beretrieved, in the order in which the fields are to be returned

        :returns: 2D array of data to be returned which can be accessed by data[field_index][datapoint]

        """
        data = []
        for field in fields:
            if isinstance(field, int):
                data.append(array(self.data.values())[field])
            elif field in self.data.keys():
                data.append(array(self.data[field]))
            else:
                data.append(None)
        return data


# often fd/4 is used in literature: check
# this is fd or cf: darcy friction factor
def fd(Re):
    if Re < 2000:
        return 64 / Re
    elif Re > 4000:
        return 1 / real(0.838 * lambertw(0.629 * Re)) ** 2
    else:
        laminar = 64 / Re
        turbulent_smooth = 1 / real(0.838 * lambertw(0.629 * Re)) ** 2
        return (laminar * (4000 - Re) + turbulent_smooth * (Re - 2000)) / (4000 - 2000)


def dp_annular_gap(D_outer, D_inner, mdot, L, rho, mu, dp=0.0):
    """Calculates the frictional and accelerationsal pressure drop over an annular gap"""
    A = np.pi * (D_outer ** 2 - D_inner ** 2) / 4
    Dh = D_outer - D_inner
    V = mdot / (rho * A)

    Re = rho * V * Dh / mu
    return 0.5 * rho * (1 + fd(Re) * L / Dh) * V ** 2 - dp


def dp_ipa_throttle(A_rat, D_outer, D_inner, mdot, rho, dp=0.0, cd=1.0):
    A = A_rat * np.pi * (D_outer ** 2 - D_inner ** 2) / 4
    V = mdot / (rho * A)
    return 0.5 * rho * V ** 2 * (1 - A_rat ** 2) / cd ** 2 - dp
