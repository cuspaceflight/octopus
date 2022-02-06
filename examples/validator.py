import os
import sys
import argparse
from pathlib import Path

import h5py
import numpy as np
from matplotlib import pyplot as plt


def main():
    # I don't think any of this is needed as it's in the HD5 already
    dirname = 'G:\\Shared drives\\Cambridge University Spaceflight\\rclone is the shit and so is robin\\Rockets\\Hybrid\\Static Firing 01 19'

    data = {}
    f = h5py.File('testfire.hd5', 'a')
    filenames = []
    for filename in os.listdir(dirname):
        if filename.split('.')[-1] == 'csv' and filename[:2] != 'PT':
            filenames.append(filename.split('.')[0])
    print(filenames)
    for filename in np.sort(filenames):
        pathname = os.path.join(dirname, filename + '.csv')
        data[filename] = np.loadtxt(pathname, delimiter=',')
        print(f'loaded {filename}')
        f.create_dataset(filename, data=data[filename])
        print(f'saved {filename} to hd5')


def load_file(path, type):

    if type == "HD5":
        hd5_file = h5py.File(path, 'r')
        print("HD5 file loaded")
        print("Found following channels")
        for channel in hd5_file['channels']: print(channel)

        #print(list(hd5_file['channels']))

        res = []
        grp = 'channels'
        for dset in hd5_file[grp]:
            if dset[:2] == 'PT':
                res.append(hd5_file[f'{grp}/{dset}'])
        return res
    else:
        raise ValueError(f"Unknown file type: {type}")


def plot_datasets(datasets, max_datapoints):

    plt.figure(0)
    for d in datasets:
        # Reduce data to not kill matplotlib
        # Ideally this is either set dynamically for a zoom level
        # Or we use some other plotting library that does not suck

        stepsize = 1 if len(d['time']) < max_datapoints else len(d['time']) // max_datapoints
        plt.plot(d['time'][::stepsize], d['data'][::stepsize], label=d)

    plt.legend(loc='upper left')
    plt.show()


if __name__ == "__main__":

    # Process arguments
    parser = argparse.ArgumentParser("Compares test fire data with octopus simulations")
    parser.add_argument("--file", "-f", type=Path, help="File to load data from")
    parser.add_argument("--type", type=str, help="Specify filetype: HD5 [more to be added??]", default="HD5")
    parser.add_argument("--max_datapoints", type=int, help="Max datapoints for plotting", default=5000)

    args = parser.parse_args()

    if args.file is None:
        print(f"File not specified, exiting...", file=sys.stderr)
        exit(1)

    # args.type should never be none due to the default value
    # if it is, then the python standard library is broken anyway
    datasets = load_file(args.file, args.type)

    plot_datasets(datasets, args.max_datapoints)

    #openhd5()
