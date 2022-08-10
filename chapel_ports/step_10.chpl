use util;

// define default simulation parameters
config const nx = 50,
             ny = 50,
             num_iterations = 100;

const dx = 2.0 / (nx - 1),
      dy = 1.0 / (ny - 1);

// create 2D arrays to represent the solution and source
const cdom = {0..<nx, 0..<ny};
const cdom_inner : subdomain(cdom) = cdom.expand((-1, -1));

var p : [cdom] real,
    b : [cdom] real;

// set source
b = 0.0;
b[nx/4, ny/4] = 100.0;
b[3*nx/4, 3*ny/4] = -100.0;

// solve the 2D Laplace's equation with the given boundary conditions
solvePoisson2D(p, b, new zeroBoundary(), dx, dy, num_iterations);

// write data
write_array_to_file("./sim_output/step_10/ch_u.txt", p);
write_array_to_file("./sim_output/step_10/ch_x.txt", linspace(0.0, 2.0, nx));
write_array_to_file("./sim_output/step_10/ch_y.txt", linspace(0.0, 1.0, ny));

// procedure to solve Laplace's Equation on p with the desired tolerance
proc solvePoisson2D(
    ref p : [?Dp] real,
    ref b : [?Db] real,
    boundaryCond,
    dx, dy,
    num_iters
) where Dp.rank == 2 && Db.rank == 2 {
    // initialize p if not already
    p = 0.0;
    boundaryCond(p);

    // create a temporary copy of p
    var pn = p;

    // iteratively solve for num_iters iterations
    for it in 0..#num_iters {
        p <=> pn;

        foreach (i, j) in cdom {
            p[i, j] = (
                dx**2 * (pn[i+1, j] + pn[i-1, j]) +
                dy**2 * (pn[i, j+1] + pn[i, j-1]) -
                b[i, j] * dx**2 * dy**2
            ) / (2.0 * (dx**2 + dy**2));
        }

        // apply boundary conditions
        boundaryCond(p);
    }
}

// procedure to apply the given boundary conditions to p
record zeroBoundary {
    proc this(ref p : [?D] real) {
        p[.., 0] = 0.0;
        p[.., D.dim(1).high] = 0.0;
        p[0, ..] = 0.0;
        p[D.dim(0).high, ..] = 0.0;
    }
}
