use StencilDist;
use util;

// define default simulation parameters
config const nx = 41,
             ny = 41,
             nt = 120,
             nu = 0.01,
             sigma = 0.0009;

const dx = 2.0 / (nx - 1),
      dy = 2.0 / (ny - 1),
      dt = sigma * dx * dy / nu;

config const write_data = false;

writeln("Running 2D Burgers' Eq. Simulation over: ");
writeln();
writeln("0 \t\t ", dx * (nx - 1));
writeln("*-----(", nx, "x)-----* \t (dx = ", dx, ")");
writeln("|\t\t|");
writeln("|\t\t|");
writeln(ny, "y \t (dy = ", dy, ")");
writeln("|\t\t|");
writeln("|\t\t|");
writeln("*---------------*");
writeln(dy * (ny - 1));
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");
writeln("with: nu = ", nu);

// setup a stencil-optimized 2D domain map for an efficient memory-parallel computation
const cdom = {0..<nx, 0..<ny};
const CDOM = cdom dmapped Stencil(cdom.expand((-1, -1)), fluff=(1, 1));
const CDOM_INNER: subdomain(CDOM) = CDOM.expand((-1, -1));

// create two 2-dimensional arrays to represent the solution in each direction
var u : [CDOM] real;
var v : [CDOM] real;

// set up the initial conditions
u = 1.0;
u[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;
u.updateFluff();

v = 1.0;
v[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;
v.updateFluff();

// apply the fd equation for nt iterations
var un : [CDOM] real = u;
var vn : [CDOM] real = v;
for i in 0..#nt {
    u <=> un;
    v <=> vn;

    // update the cached "fluff" points on the edge of each locale
    un.updateFluff();
    vn.updateFluff();

    forall (i, j) in CDOM_INNER {
        u[i, j] = un[i, j]
                - (dt / dx * un[i, j] * (un[i, j] - un[i-1, j]))
                - (dt / dy * vn[i, j] * (un[i, j] - un[i, j-1]))
                + (nu * dt / dx**2 * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]))
                + (nu * dt / dy**2 * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1]));


        v[i, j] = vn[i, j]
                - (dt / dx * un[i, j] * (vn[i, j] - vn[i-1, j]))
                - (dt / dy * vn[i, j] * (vn[i, j] - vn[i, j-1]))
                + (nu * dt / dx**2 * (vn[i+1, j] - 2.0 * vn[i, j] + vn[i-1, j]))
                + (nu * dt / dy**2 * (vn[i, j+1] - 2.0 * vn[i, j] + vn[i, j-1]));
    }

    u[0, ..] = 1.0;
    u[.., 0] = 1.0;
    u[nx - 1, ..] = 1.0;
    u[.., ny - 1] = 1.0;

    v[0, ..] = 1.0;
    v[.., 0] = 1.0;
    v[nx - 1, ..] = 1.0;
    v[.., ny - 1] = 1.0;
}

if write_data {
    write_array_to_file("./sim_output/step_8/ch_u.txt", u);
    write_array_to_file("./sim_output/step_8/ch_v.txt", v);
    write_array_to_file("./sim_output/step_8/ch_x.txt", linspace(0.0, 2.0, nx));
    write_array_to_file("./sim_output/step_8/ch_y.txt", linspace(0.0, 2.0, ny));
}
