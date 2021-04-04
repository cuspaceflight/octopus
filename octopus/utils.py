from collections import OrderedDict
from csv import reader

from numpy import array, nan_to_num


def derivative(f, axis, *args, dx=0.01):
    """Calculate the  numerical derivative of a scalar funtion f(*args) with respect to args[axis].

        Parameters:
            f (Callable): vector function to find derivative of with the signature f(*args) -> iterable
            axis (int): index of the parameter against which the derivative of f is to be taken
            args (float): parameters to pass to f, including the one against which the derivative is to be taken
            dx (float): differentiation step
        Returns:
            derivative (float): array of floats representing vector derivative f with respect to args[axis]
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

    def __init__(self, filename):
        """Retrieve thermodynamic data from a tab-delimited csv file

        Parameters:
            filename (str): filename, not including file extension, of a tab-delimited NIST datafile"""
        with open(f'{filename}.csv', newline='') as f:
            csvreader = reader(f, delimiter='\t')
            rows = [row for row in csvreader]
            self.data = OrderedDict()
            for i, header in enumerate(rows[0]):
                self.data[header] = [float(row[i]) if row[i] != 'undefined' else None for row in rows[1:]]

    def list_fields(self):
        """List fields in NIST datafile

        Returns:
            keys (str): list of column headers in the NIST datafile"""
        return [key for key in self.data.keys()]

    def get_fields(self, *fields):
        """Retrieve data under specified fields

        Parameters:
            fields (str or int): titles or indices of column headers in the datafile for which data is to be
            retrieved, in the order in which the fields are to be returned

        Returns:
            data (float): 2D array of data to be returned which can be accessed by data[field_index][datapoint]"""
        data = []
        for field in fields:
            if isinstance(field, int):
                data.append(array(self.data.values())[field])
            elif field in self.data.keys():
                data.append(array(self.data[field]))
            else:
                data.append(None)
        return data
