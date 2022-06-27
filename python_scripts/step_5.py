import numpy
import sympy
from matplotlib import pyplot
from matplotlib import cm
import sys

show_plots = "--show_plots" in sys.argv[1:]

### setup simulation parameters
nx = 81
ny = 81
nt = 100
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .2
dt = sigma * dx

### setup computational domain
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))

### setup initial conditions
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2

fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
if show_plots: pyplot.show()

### apply fd equation for nt iterations
for n in range(nt):
    un = u.copy()
    u[1:, 1:] = (un[1:, 1:] - (c * dt / dx * (un[1:, 1:] - un[1:, :-1])) -
                              (c * dt / dy * (un[1:, 1:] - un[:-1, 1:])))
    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
if show_plots: pyplot.show()

numpy.savetxt("./sim_output/step_5_py_output.txt", u, fmt='%.8f')
