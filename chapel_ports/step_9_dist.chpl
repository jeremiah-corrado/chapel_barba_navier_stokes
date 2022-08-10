use StencilDist;
use BlockDist;
use util;

// define default simulation parameters
config const nx = 31,
             ny = 31,
             l1_tolerance = 1e-4,
             max_num_iters = 10000;

const dx = 2.0 / (nx - 1),
      dy = 1.0 / (ny - 1);

config const write_data = false;

const cdom = {0..<nx, 0..<ny};
const CDOM = cdom dmapped Stencil(cdom.expand((-1, -1)), fluff=(1, 1));
const CDOM_INNER: subdomain(CDOM) = CDOM.expand((-1, -1));

// create a 2D array to represent solution
var p: [CDOM] real;

// solve the 2D Laplace's equation with the given boundary conditions
solveLaplace2D(p, new linYBoundary(), dx, dy, l1_tolerance);

if write_data {
    write_array_to_file("./sim_output/step_9/ch_u.txt", p);
    write_array_to_file("./sim_output/step_9/ch_x.txt", linspace(0.0, 2.0, nx));
    write_array_to_file("./sim_output/step_9/ch_y.txt", linspace(0.0, 1.0, ny));
}

// procedure to solve Laplace's Equation on p with the desired tolerance
proc solveLaplace2D(
    ref p : [?D] real,
    boundaryCond,
    dx, dy,
    l1_rel_delta_tolerance
) where D.rank == 2 {
    // initialize p if not already
    p = 0.0;
    boundaryCond(p);

    // define some working variables
    var l1norm_rel_delta = 1.0;
    var pn : [CDOM] real = p;
    var i = 0;

    // iteratively solve until desired tolerance is reached
    while l1norm_rel_delta > l1_rel_delta_tolerance && i < max_num_iters {
        p <=> pn;

        // apply fd equation
        forall (i, j) in CDOM_INNER {
            p[i, j] = (
                dy**2 * (pn[i, j+1] + pn[i, j-1]) +
                dx**2 * (pn[i+1, j] + pn[i-1, j])
            ) / (2.0 * (dx**2 + dy**2));
        }

        // apply boundary conditions
        boundaryCond(p);

        // compute delta in L1 norm
        l1norm_rel_delta = (+ reduce (abs(p) - abs(pn))) / (+ reduce abs(pn));

        i +=1;
    }

    writeln("Ran for ", i, " iterations");
}

// procedure to apply the given boundary conditions to p
record linYBoundary {
    const cdom_y = {0..<ny};
    const CDOM_Y = cdom_y dmapped Block(boundingBox = cdom_y.expand(-1));
    var y = linspace_dist(0, 1.0, ny, CDOM_Y);

    proc this(ref p : [?D] real) {
        p[.., 0] = 0.0;                                     // p(0.0, y) = 0
        p[.., D.dim(1).high] = this.y;                      // p(2.0, y) = y
        p[0, ..] = p[1, ..];                                // dp/dy(x, 0.0) = 0
        p[D.dim(0).high, ..] = p[D.dim(0).high-1, ..];      // dp/dy(x, 1.0) = 0

        p.updateFluff();
    }
}
