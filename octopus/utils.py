from collections import OrderedDict
from csv import reader

from numpy import array, nan_to_num


def derivative(f, axis, *args, dx=0.01):
    if axis >= len(args):
        raise ValueError("axis index too high")
    args_l = list(args)
    args_h = list(args)
    args_l[axis] -= dx
    args_h[axis] += dx

    return nan_to_num(array((f(*args_h) - f(*args_l)) / (2 * dx)))


class Nist:
    def __init__(self, filename):
        with open(f'{filename}.csv', newline='') as f:
            csvreader = reader(f, delimiter='\t')
            rows = [row for row in csvreader]
            self.data = OrderedDict()
            for i, header in enumerate(rows[0]):
                self.data[header] = [float(row[i]) if row[i] != 'undefined' else None for row in rows[1:]]

    def list_fields(self):
        return [key for key in self.data.keys()]

    def get_fields(self, *fields):
        data = []
        for field in fields:
            if isinstance(field, int):
                data.append(array(self.data.values())[field])
            elif field in self.data.keys():
                data.append(array(self.data[field]))
            else:
                data.append(None)
        return data
