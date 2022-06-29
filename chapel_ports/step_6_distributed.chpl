use StencilDist;
use util;

// define default simulation parameters
config const nx = 81;
config const ny = 81;
config const nt = 100;
config const c = 1;
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
config const sigma = 0.2;
const dt = sigma * dx;

writeln("Running 2D Non-Linear Convection Simulation over: ");
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
writeln("with: c = ", c);

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
v = 1.0;
v[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;

// apply the fd equation for nt iterations
var un = u;
var vn = v;
for i in 0..#nt {
    u <=> un;
    v <=> vn;

    // compute the stencil computation in parallel across all locales
    forall (i, j) in SpaceInner {
        u[i, j] = un[i, j] - (un[i, j] * c * dt / dx * (un[i, j] - un[i, j - 1])) -
                             (vn[i, j] * c * dt / dy * (un[i, j] - un[i - 1, j]));

        v[i, j] = vn[i, j] - (un[i, j] * c * dt / dx * (vn[i, j] - vn[i, j - 1])) -
                             (vn[i, j] * c * dt / dy * (vn[i, j] - vn[i - 1, j]));
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

write_array_to_file("./sim_output/step_6_dist_u_output.txt", u);
write_array_to_file("./sim_output/step_6_dist_v_output.txt", v);
