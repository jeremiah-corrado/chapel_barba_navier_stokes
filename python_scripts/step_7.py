import numpy
from matplotlib import pyplot
from matplotlib import cm
import sys
import math

show_plots = "--show_plots" in sys.argv[1:]

### setup simulation parameters
nx = 31
ny = 31
nt = 17
nu = .05
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .25
dt = sigma * dx * dy / nu

### setup computational domain
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
un = numpy.ones((ny, nx))

### define diffusion function
def diffuse(nt):
    ## set initial conditions
    u[:] = 1
    u[int(.5 / dy):int(1 / dy + 1), int(.5 / dx):int(1 / dx + 1)] = 2

    ## apply fd equation
    for n in range(nt):
        un = u.copy()

        u[1:-1, 1:-1] = (un[1:-1,1:-1] +
                            nu * dt / dx**2 *
                            (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                            nu * dt / dy**2 *
                            (un[2:,1: -1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

        u[0, :] = 1
        u[-1, :] = 1
        u[:, 0] = 1
        u[:, -1] = 1

    ## generate plot
    fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
    X, Y = numpy.meshgrid(x, y)
    surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
    ax.set_title("u(x, y)")
    if show_plots: pyplot.show()

diffuse(10)
numpy.savetxt("./sim_output/step_7_a_py_output.txt", u, fmt='%.8f')

diffuse(14)
numpy.savetxt("./sim_output/step_7_b_py_output.txt", u, fmt='%.8f')

diffuse(50)
numpy.savetxt("./sim_output/step_7_py_output.txt", u, fmt='%.8f')
