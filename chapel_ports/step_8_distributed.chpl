use StencilDist;
use util;

// define default simulation parameters
config const nx = 41;
config const ny = 41;
config const nt = 120;
config const nu = 0.01;
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
config const sigma = 0.0009;
const dt = sigma * dx * dy / nu;

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
const Space = {0..<nx, 0..<ny};
const SpaceInner = {1..<(nx-1), 1..<(nx-1)};

const CompDom = Space dmapped Stencil(
    SpaceInner, // our stencil computation is concerned with the inner set of points
    fluff=(1,1) // each locale only needs to know about 1 point from the adjacent locales
);

// create two 2-dimensional arrays to represent the solution in each direction
var u : [CompDom] real;
var v : [CompDom] real;

// set up the initial conditions
u = 1.0;
u[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;
u.updateFluff();

v = 1.0;
v[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;
v.updateFluff();

// apply the fd equation for nt iterations
var un = u;
var vn = v;
for i in 0..#nt {
    u <=> un;
    v <=> vn;

    // update the cached "fluff" points on the edge of each locale
    un.updateFluff();
    vn.updateFluff();

    forall (i, j) in SpaceInner {
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
