import numpy
import sys

show_plots = "--show_plots" in sys.argv[1:]

nx = 41
dx = 2 / (nx-1)
nt = 25
dt = .025

u = numpy.ones(nx)
u[int(.5 / dx):int(1 / dx + 1)] = 2
print("u(t = 0): \t", u);

if show_plots: 
    from matplotlib import pyplot
    pyplot.plot(numpy.linspace(0, 2, nx), u)
    pyplot.show()

un = numpy.ones(nx)

for n in range(nt):
    un = u.copy()
    for i in range(1, nx):
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1])

if show_plots: 
    pyplot.plot(numpy.linspace(0, 2, nx), u)
    pyplot.show()

print("u(t = ", nt * dt, "): \t", u)

numpy.savetxt("./sim_output/step_2/py_u.txt", u, fmt='%.8f')
