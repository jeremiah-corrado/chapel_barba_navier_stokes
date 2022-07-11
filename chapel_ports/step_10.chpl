use util;

// define default simulation parameters
config const nx = 50;
config const ny = 50;
const dx = 2.0 / (nx - 1);
const dy = 1.0 / (ny - 1);
config const num_iterations = 100;

// create a 2D array to represent solution
var p : [{0..<nx, 0..<ny}] real;

// solve the 2D Laplace's equation with the given boundary conditions
solveLaplace2D(p, new zeroBoundary(), dx, dy, num_iterations);
write_array_to_file("./sim_output/step_9/ch_u.txt", p);
write_array_to_file("./sim_output/step_9/ch_x.txt", linspace(0.0, 2.0, nx));
write_array_to_file("./sim_output/step_9/ch_y.txt", linspace(0.0, 2.0, ny));

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

    var pn = p;

    // define a subdomain on the interior points of the domain
    var Ds : subdomain(Dp); //this is not strictly necessary; however it does reduce the requisite number of bound checks during the loop below
    Ds = D[D.dim(0).expand(-1), D.dim(1).expand(-1)];

    // iteratively solve until desired tolerance is reached
    for i in 0..#num_iters {
        p <=> pn;

        // // apply fd equation
        // foreach (i, j) in Ds {
        //     p[i, j] = (
        //         dy**2 * (pn[i, j+1] + pn[i, j-1]) +
        //         dx**2 * (pn[i+1, j] + pn[i-1, j])
        //     ) / (2.0 * (dx**2 + dy**2));
        // }

        // apply boundary conditions
        boundaryCond(p);
    }

    writeln("Ran for ", i, " iterations");
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
