import numpy
from matplotlib import pyplot
import sys

show_plots = "--show_plots" in sys.argv[1:]

nx = 41
dx = 2 / (nx-1)
nt = 20
nu = 0.3
sigma = 0.2
dt = sigma * dx**2 / nu

u = numpy.ones(nx)
u[int(.5 / dx):int(1 / dx + 1)] = 2
print("u(t = 0): \t", u);

pyplot.plot(numpy.linspace(0, 2, nx), u)
if show_plots: pyplot.show()

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx - 1):
        u[i] = un[i] + nu * dt / dx**2 * (un[i+1] - 2 * un[i] + un[i-1])

pyplot.plot(numpy.linspace(0, 2, nx), u)
if show_plots: pyplot.show()

print("u(t = ", nt * dt, "): \t", u)

numpy.savetxt("./sim_output/step_3_py_output.txt", u, fmt='%.8f')
