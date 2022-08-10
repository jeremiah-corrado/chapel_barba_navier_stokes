use StencilDist;
use util;

// define u-function directly
var ufunc = lambda(t:real, x:real, nu:real) {
    return -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) -
        (-8*t + 2*x - 4*pi)*exp(-(-4*t + x - 2*pi)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 2*pi)**2/(4*nu*(t + 1))) +
        exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4;
};

// define default simulation parameters
config const nx = 101;
config const nt = 100;
const dx = 2 * pi / (nx - 1);
config const nu = 0.07;
const dt = dx * nu;

config const write_data = false;

writeln("Running 1D Diffusion Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t  ", dx * (nx - 1));
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");

// setup a stencil-optimized domain map for an efficient memory-parallel computation
const cdom = 0..<nx;
const CDOM = cdom dmapped Stencil(cdom.expand(-1), fluff=(1,));
const CDOM_INNER : subdomain(CDOM) = cdom.expand(-1);

var u : [CDOM] real;

// setup initial conditions
const x = linspace_dist(0.0, 2 * pi, nx, CDOM);
[i in Space] u[i] = ufunc(0, x[i], nu);
u.updateFluff();

// apply the fd equation for nt iterations
var un : [CDOM] real = u;
for n in 0..#nt {
    u <=> un;
    un.updateFluff(); // update the cached "fluff" points on the edge of each locale

    // compute the bulk of the stencil computation in parallel across all locales
    forall i in CDOM_INNER {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1]) + nu * dt / dx**2 *
                (un[i+1] - 2 * un[i] + un[i-1]);
    }

    // apply the cyclic boundary condition
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx-2]) + nu * dt / dx**2 *
                (un[1] - 2 * un[0] + un[nx-2]);
    u[nx-1] = u[0];
}

if write_data {
    write_array_to_file("./sim_output/step_4/ch_u.txt", u);
    write_array_to_file("./sim_output/step_4/ch_x.txt", x);
}
