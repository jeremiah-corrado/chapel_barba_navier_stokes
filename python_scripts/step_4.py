import numpy
import sympy
import sys

show_plots = "--show_plots" in sys.argv[1:]

### define function for computing initial conditions
x, nu, t = sympy.symbols('x nu t')
phi = (sympy.exp(-(x - 4 * t)**2 / (4 * nu * (t + 1))) +
       sympy.exp(-(x - 4 * t - 2 * sympy.pi)**2 / (4 * nu * (t + 1))))
print("phi: \t", phi)

phiprime = phi.diff(x)
print("phi': \t", phiprime)

from sympy.utilities.lambdify import lambdify

u = -2 * nu * (phiprime / phi) + 4
print("u(t, x, nu): \t", u)

ufunc = lambdify((t, x, nu), u)
print("u(1, 4, 3): \t", ufunc(1, 4, 3))

### setup simulation parameters
nx = 101
nt = 100
dx = 2 * numpy.pi / (nx - 1)
nu = .07
dt = dx * nu

### setup initial conditions
x = numpy.linspace(0, 2 * numpy.pi, nx)
un = numpy.empty(nx)
t = 0

u = numpy.asarray([ufunc(t, x0, nu) for x0 in x])
print("u(t = 0, x): \t", u)

if show_plots: 
    from matplotlib import pyplot
    pyplot.figure(figsize=(11, 7), dpi=100)
    pyplot.plot(x, u, marker='o', lw=2)
    pyplot.xlim([0, 2 * numpy.pi])
    pyplot.ylim([0, 10])
    pyplot.show()

### Apply Differential Equation
for n in range(nt):
    un = u.copy()
    for i in range(1, nx-1):
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 *\
                (un[i+1] - 2 * un[i] + un[i-1])
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[-2]) + nu * dt / dx**2 *\
                (un[1] - 2 * un[0] + un[-2])
    u[-1] = u[0]

u_analytical = numpy.asarray([ufunc(nt * dt, xi, nu) for xi in x])


if show_plots: 
    pyplot.figure(figsize=(11, 7), dpi=100)
    pyplot.plot(x,u, marker='o', lw=2, label='Computational')
    pyplot.plot(x, u_analytical, label='Analytical')
    pyplot.xlim([0, 2 * numpy.pi])
    pyplot.ylim([0, 10])
    pyplot.legend()
    pyplot.show()

numpy.savetxt("./sim_output/step_4/py_u.txt", u, fmt='%.8f')
