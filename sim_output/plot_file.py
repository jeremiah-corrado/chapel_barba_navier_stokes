from matplotlib import pyplot
from matplotlib import cm
import numpy
import sys

data_file_name = sys.argv[1]
print("Plotting: ", data_file_name)

plot_data = numpy.loadtxt(data_file_name)

if plot_data.ndim == 1:
    x = numpy.linspace(0, 2.0, plot_data.shape[0])
    pyplot.plot(x, plot_data)

elif plot_data.ndim == 2:
    x = numpy.linspace(0, 2.0, plot_data.shape[0])
    y = numpy.linspace(0, 2.0, plot_data.shape[1])

    fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
    X, Y = numpy.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, plot_data[:], cmap=cm.viridis)
else:
    print("no")

pyplot.show()
