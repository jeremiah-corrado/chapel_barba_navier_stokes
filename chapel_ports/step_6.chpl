use util;

// define default simulation parameters
config const nx = 81,
             ny = 81,
             nt = 100,
             c = 1,
             sigma = 0.2;

const dx = 2.0 / (nx - 1),
      dy = 2.0 / (ny - 1),
      dt = sigma * dx;

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

// create two 2-dimensional arrays to represent the solution in each direction
var cdom = {0..<nx, 0..<ny};
var cdom_inner : subdomain(cdom) = cdom.expand((-1, -1));
var u : [cdom] real;
var v : [cdom] real;

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

    forall (i, j) in cdom_inner {
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

write_array_to_file("./sim_output/step_6/ch_u.txt", u);
write_array_to_file("./sim_output/step_6/ch_v.txt", v);
write_array_to_file("./sim_output/step_6/ch_x.txt", linspace(0.0, 2.0, nx));
write_array_to_file("./sim_output/step_6/ch_y.txt", linspace(0.0, 2.0, ny));
