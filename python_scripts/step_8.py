import numpy
import sympy
from matplotlib import pyplot
from matplotlib import cm
import sys

show_plots = "--show_plots" in sys.argv[1:]

### setup simulation parameters
nx = 41
ny = 41
nt = 120
c = 1
dx = 2 / (nx - 1)
dy = 2 / (ny - 1)
sigma = .0009
nu = 0.01
dt = sigma * dx * dy / nu

### setup computational domain
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 2, ny)

u = numpy.ones((ny, nx))
v = numpy.ones((ny, nx))

un = numpy.ones((ny, nx))
vn = numpy.ones((ny, nx))

### setup initial conditions
u[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2
v[int(.5 / dy):int(1 / dy + 1),int(.5 / dx):int(1 / dx + 1)] = 2


fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
X, Y = numpy.meshgrid(x, y)
surf = ax.plot_surface(X, Y, u[:], cmap=cm.viridis)
if show_plots: pyplot.show()

### apply fd equation for nt iterations
for n in range(nt):
    un = u.copy()
    vn = v.copy()

    u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                     dt / dx * un[1:-1, 1:-1] *
                     (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                     dt / dy * vn[1:-1, 1:-1] *
                     (un[1:-1, 1:-1] - un[0:-2, 1:-1]) +
                     nu * dt / dx**2 *
                     (un[1:-1,2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                     nu * dt / dy**2 *
                     (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1]))

    v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                     dt / dx * un[1:-1, 1:-1] *
                     (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                     dt / dy * vn[1:-1, 1:-1] *
                    (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) +
                     nu * dt / dx**2 *
                     (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                     nu * dt / dy**2 *
                     (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1]))

    u[0, :] = 1
    u[-1, :] = 1
    u[:, 0] = 1
    u[:, -1] = 1

    v[0, :] = 1
    v[-1, :] = 1
    v[:, 0] = 1
    v[:, -1] = 1

fig, ax = pyplot.subplots(1, 2, subplot_kw={"projection": "3d"})
X, Y = numpy.meshgrid(x, y)
surf = ax[0].plot_surface(X, Y, u[:], cmap=cm.viridis)
ax[0].set_title("u(x, y)")
surf = ax[1].plot_surface(X, Y, v[:], cmap=cm.viridis)
ax[1].set_title("v(x, y)")
if show_plots: pyplot.show()

numpy.savetxt("./sim_output/step_8/py_u.txt", u, fmt='%.8f')
numpy.savetxt("./sim_output/step_8/py_v.txt", v, fmt='%.8f')
