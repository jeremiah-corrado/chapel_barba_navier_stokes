
// define a default setup
config const nx = 41;
const dx : real = 2.0 * pi / (nx - 1);
config const nt = 25;
config const nu = 0.3;
config const sigma = 0.2;
const dt = sigma * dx**2 / nu;

writeln("Running 1D Burgers Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t  ", dx * nx);
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");
writeln("with: nu = ", nu);

// create an array to represent the computational Domain
const Space = {0..nx};
const SpaceInner = {1..nx-1};
var u : [Space] real;

// set up the initial conditions
u = 1.0;
u[(0.5 / dx):int..(1.0 / dx + 1):int] = 2.0;

writeln("\nDomain (t = 0):");
writeln(u);

// apply the differential equation for nt iterations
var un : [Space] real;

for n in 0..nt {
    u <=> un;

    foreach i in SpaceInner {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / dx**2 *
            (un[i + 1] - 2 * un[i] + un[i - 1]);
    }

    // cyclical boundary conditions
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 2]) + nu * dt / dx**2 *
        (un[1] -  2 * un[0] + un[nx - 2]);
    u[nx - 1] = u[0];
}

writeln("\nDomain (t = ", nt * dt,"):");
writeln(u);
