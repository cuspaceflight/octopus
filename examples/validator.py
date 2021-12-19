import os

import h5py
import numpy as np
from matplotlib import pyplot as plt


def main():
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


def openhd5():
    f = h5py.File('20190122-006.h5', 'r')
    print(list(f['channels']))

    plt.figure(0)
    grp = 'channels'
    for dset in f[grp]:
        if dset[:2]=='TC':
            d = f[f'{grp}/{dset}']
            plt.plot(d['time'][:], d['data'][:], label=dset)
    plt.legend(loc='upper left')
    plt.show()


if __name__ == "__main__":
    openhd5()
