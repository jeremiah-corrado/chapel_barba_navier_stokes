import numpy
import sympy
from matplotlib import pyplot
from matplotlib import cm
import sys

# Parameters
nx = 50
ny = 50
nt  = 100
xmin = 0
xmax = 2
ymin = 0
ymax = 1

dx = (xmax - xmin) / (nx - 1)
dy = (ymax - ymin) / (ny - 1)

show_plots = "--show_plots" in sys.argv[1:]

# Initialization
p  = numpy.zeros((ny, nx))
pd = numpy.zeros((ny, nx))
b  = numpy.zeros((ny, nx))
x  = numpy.linspace(xmin, xmax, nx)
y  = numpy.linspace(ymin, ymax, ny)

# Source
b[int(ny / 4), int(nx / 4)]  = 100
b[int(3 * ny / 4), int(3 * nx / 4)] = -100

# Iteratively Solve Poisson's Equation
for it in range(nt):
    pd = p.copy()

    p[1:-1,1:-1] = (((pd[1:-1, 2:] + pd[1:-1, :-2]) * dy**2 +
                    (pd[2:, 1:-1] + pd[:-2, 1:-1]) * dx**2 -
                    b[1:-1, 1:-1] * dx**2 * dy**2) /
                    (2 * (dx**2 + dy**2)))

    p[0, :] = 0
    p[ny-1, :] = 0
    p[:, 0] = 0
    p[:, nx-1] = 0

### define a plotting function
def plot2D(x, y, p):
    if show_plots:
        fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
        X, Y = numpy.meshgrid(x, y)
        surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
                linewidth=0, antialiased=False)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)
        ax.view_init(30, 225)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        pyplot.show();


plot2D(x, y, p)
