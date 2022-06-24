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

writeln("Running 1D Diffusion Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t  ", dx * (nx - 1));
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");

// setup a stencil-optimized domain map for an efficient memory-parallel computation
const Space = {0..<nx};
const SpaceInner = {1..<(nx-1)};

const CompDom = Space dmapped Stencil(
    SpaceInner, // our stencil computation is concerned with the inner set of points
    fluff=(1,) // each local only needs to know about 1 point from the adjacent locals
);
var u : [CompDom] real;

// setup initial conditions
const x = linspace(0.0, 2 * pi, nx);
[i in Space] u[i] = ufunc(0, x[i], nu);
u.updateFluff();

writeln("u(t = 0, x):");
writeln(u);

// apply the differential equation for nt iterations
var un = u;
for n in 0..#nt {
    u <=> un;
    un.updateFluff(); // update the cached "fluff" points on the edge of each local

    // compute the bulk of the stencil computation in parallel across all locals
    forall i in CompDom {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1]) + nu * dt / dx**2 *
                (un[i+1] - 2 * un[i] + un[i-1]);
    }

    // apply the cyclic boundary condition
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx-2]) + nu * dt / dx**2 *
                (un[1] - 2 * un[0] + un[nx-2]);
    u[nx-1] = u[0];
}

writeln("Domain (t = ", nt * dt,"):");
writeln(u[Space]);

write_array_to_file("./sim_output/step_4_dist_output.txt", u);
