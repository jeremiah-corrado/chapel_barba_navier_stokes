from matplotlib import pyplot
from matplotlib import cm
import numpy
import sys
import os
import re

folder_name = sys.argv[1]
print("Plotting From: ", folder_name)

ch_files = dict()
for f in os.listdir(folder_name):
    x = re.search("^ch_(.*)\.txt$", f)
    if x:
        ch_files[x.group(1)] = numpy.loadtxt(folder_name + "/" + f)

if not 'u' in ch_files or not 'x' in ch_files:
    print("ensure data files are in given folder: ", folder_name);
else:
    if ch_files['u'].ndim == 1:
        fig, ax = pyplot.subplots()
        ax.set_xlabel('$x$')
        ax.set_title('u(x)')
        pyplot.plot(ch_files['x'], ch_files['u'])

    elif ch_files['u'].ndim == 2:
        X, Y = numpy.meshgrid(ch_files['x'], ch_files['y'])
        if 'v' in ch_files:
            fig, ax = pyplot.subplots(1, 2, subplot_kw={"projection": "3d"})
            surf = ax[0].plot_surface(X, Y, ch_files['u'], cmap=cm.viridis)
            ax[0].set_title('u(x, y)')
            ax[0].set_xlabel('$x$')
            ax[0].set_ylabel('$y$')
            surf = ax[1].plot_surface(X, Y, ch_files['v'], cmap=cm.viridis)
            ax[1].set_title('v(x, y)')
            ax[1].set_xlabel('$x$')
            ax[1].set_ylabel('$y$')
        else:
            fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
            surf = ax.plot_surface(X, Y, ch_files['u'][:], cmap=cm.viridis)
            ax.set_title('u(x, y)')
            ax.set_xlabel('$x$')
            ax.set_ylabel('$y$')
    else:
        print("Only 1 or 2 dimension plotting is supported")

    pyplot.show()
