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

if 'p' in ch_files and 'u' in ch_files and 'v' in ch_files and 'x' in ch_files and 'y' in ch_files:
    X, Y = numpy.meshgrid(ch_files['x'], ch_files['y'])
    p = ch_files['p']

    fig = pyplot.figure(figsize=(11,7), dpi=100)
    pyplot.contourf(X, Y, p, alpha=0.5, cmap=cm.viridis)
    pyplot.colorbar()

    # plotting the pressure field outlines
    pyplot.contour(X, Y, p, cmap=cm.viridis)
    pyplot.streamplot(X, Y, ch_files['u'], ch_files['v'])
    pyplot.xlabel('X')
    pyplot.ylabel('Y')
elif 'u' in ch_files and 'x' in ch_files:
    if ch_files['u'].ndim == 1:
        fig, ax = pyplot.subplots()
        ax.set_xlabel('$x$')
        ax.set_title('u(x)')
        pyplot.plot(ch_files['x'], ch_files['u'])

    elif ch_files['u'].ndim == 2 and 'y' in ch_files:
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
else:
    print("ensure all data files are in given folder: ", folder_name);

pyplot.show()
