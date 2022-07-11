import numpy
import sympy
from matplotlib import pyplot
from matplotlib import cm
import sys

show_plots = "--show_plots" in sys.argv[1:]

### setup simulation parameters
nx = 31
ny = 31
c = 1
dx = 2 / (nx - 1)
dy = 1 / (ny - 1)

### define a function to solve Laplace's Equation
def laplace2d(p, y, dx, dy, l1norm_target):
    l1norm = 1
    pn = numpy.empty_like(p)

    i = 0

    while l1norm > l1norm_target:
        pn = p.copy()
        p[1:-1, 1:-1] = ((dy**2 * (pn[1:-1, 2:] + pn[1:-1, 0:-2]) +
                         dx**2 * (pn[2:, 1:-1] + pn[0:-2, 1:-1])) /
                        (2 * (dx**2 + dy**2)))

        p[:, 0] = 0  # p = 0 @ x = 0
        p[:, -1] = y  # p = y @ x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
        p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

        print(numpy.sum(numpy.abs(p)));

        l1norm = (
            numpy.sum(numpy.abs(p[:]) - numpy.abs(pn[:])) /
            numpy.sum(numpy.abs(pn[:])))

        # print(l1norm)

        i = i + 1

    print("Ran for ", i, " iterations")

    return p

### define a plotting function
def plot2D(x, y, p):
    if show_plots:
        fig, ax = pyplot.subplots(subplot_kw={"projection": "3d"})
        X, Y = numpy.meshgrid(x, y)
        surf = ax.plot_surface(X, Y, p[:], rstride=1, cstride=1, cmap=cm.viridis,
                linewidth=0, antialiased=False)
        ax.set_xlim(0, 2)
        ax.set_ylim(0, 1)
        ax.view_init(30, 225)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        pyplot.show();


### setup initial conditions
p = numpy.zeros((ny, nx))  # create a XxY vector of 0's

### setup plotting space
x = numpy.linspace(0, 2, nx)
y = numpy.linspace(0, 1, ny)

### setup boundary conditions
p[:, 0] = 0  # p = 0 @ x = 0
p[:, -1] = y  # p = y @ x = 2
p[0, :] = p[1, :]  # dp/dy = 0 @ y = 0
p[-1, :] = p[-2, :]  # dp/dy = 0 @ y = 1

print(numpy.sum(numpy.abs(p)))
print()

### Plot initial conditions
plot2D(x, y, p)

### Solve Laplace's Equation
p = laplace2d(p, y, dx, dy, 1e-4)

### Plot solution
plot2D(x, y, p)

numpy.savetxt("./sim_output/step_9/py_u.txt", p, fmt='%.8f')
