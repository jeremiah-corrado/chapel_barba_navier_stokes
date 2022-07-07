use StencilDist;
use util;

// define default simulation parameters
config const nx = 31;
config const ny = 31;
config const nu = 0.05;
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
config const sigma = 0.25;
const dt = sigma * dx * dy / nu;

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
const Space = {0..<nx, 0..<ny};
const SpaceInner = {1..<(nx-1), 1..<(nx-1)};

const CompDom = Space dmapped Stencil(
    SpaceInner, // our stencil computation is concerned with the inner set of points
    fluff=(1,1) // each locale only needs to know about 1 point from the adjacent locales
);

// create an 2-dimensional array to represent the computational Domain
var u : [CompDom] real;

diffuse(10, u); write_array_to_file("./sim_output/step_7_dist_a_output.txt", u);
diffuse(14, u); write_array_to_file("./sim_output/step_7_dist_b_output.txt", u);
diffuse(50, u); write_array_to_file("./sim_output/step_7_dist_output.txt", u);

// apply the diffusion operation to 'u' for 'nt' iterations
proc diffuse(nt: int, ref u : [?D] real) {
    // make sure 'u' is 2D
    if D.rank != 2 then halt();

    // call helper procedure to set the initial conditions
    init_conditions(u);
    u.updateFluff();

    // apply the fd equation for nt iterations
    var un = u;
    for i in 0..#nt {
        u <=> un;

        // update the cached "fluff" points on the edge of each locale
        un.updateFluff();

        forall (i, j) in SpaceInner {
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
