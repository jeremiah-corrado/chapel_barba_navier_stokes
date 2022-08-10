use StencilDist;
use util;

// define default simulation parameters
config const nx = 31,
             ny = 31,
             nu = 0.05,
             sigma = 0.25;

const dx = 2.0 / (nx - 1),
      dy = 2.0 / (ny - 1),
      dt = sigma * dx * dy / nu;

config const write_data = false;

writeln("Running 2D Diffusion Simulation over: ");
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
writeln("with dt = ", dt, ")");

// setup a stencil-optimized 2D domain map for an efficient memory-parallel computation
const cdom = {0..<nx, 0..<ny};
const CDOM = cdom dmapped Stencil(cdom.expand((-1, -1)), fluff=(1, 1));
const CDOM_INNER: subdomain(CDOM) = CDOM.expand((-1, -1));

// create an 2-dimensional array to represent the computational Domain
var u : [CDOM] real;

diffuse(50, u);

if write_data {
    write_array_to_file("./sim_output/step_7/ch_u.txt", u);
    write_array_to_file("./sim_output/step_7/ch_x.txt", linspace(0.0, 2.0, nx));
    write_array_to_file("./sim_output/step_7/ch_y.txt", linspace(0.0, 2.0, ny));
}

// apply the diffusion operation to 'u' for 'nt' iterations
proc diffuse(nt: int, ref u : [?d] real) where d.rank == 2 {

    // call helper procedure to set the initial conditions
    init_conditions(u);
    u.updateFluff();

    // apply the fd equation for nt iterations
    var un : [CDOM] real = u;
    for i in 0..#nt {
        u <=> un;

        // update the cached "fluff" points on the edge of each locale
        un.updateFluff();

        forall (i, j) in CDOM_INNER {
            u[i, j] = un[i, j] +
                        nu * dt / dx**2 *
                            (un[i-1, j] - 2 * un[i, j] + un[i+1, j]) +
                        nu * dt / dy**2 *
                            (un[i, j-1] - 2 * un[i, j] + un[i, j+1]);
        }

        u[0, ..] = 1.0;
        u[.., 0] = 1.0;
        u[nx - 1, ..] = 1.0;
        u[.., ny - 1] = 1.0;
    }
}

// set the initial condition on u
proc init_conditions(ref u : [] real) {
    u = 1.0;
    u[(0.5 / dy):int..<(1.0 / dy + 1):int, (0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;
}
