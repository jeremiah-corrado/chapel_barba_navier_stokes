
// define a default setup
config const nx = 41;
const dx : real = 2.0 / (nx - 1);
config const nt = 25;
config const dt = 0.025;

writeln("Running 1D Non-Linear Simulation with: ");
writeln("nx = ", nx, "\t dx = ", dx);
writeln("nx = ", nt, "\t dt = ", dt);
writeln();

// create an array to represent the Domain
var u : [0..#nx] real;


// set up the initial conditions
u = 1.0;
u[(0.5 / dx):int..(1.0 / dx + 1):int] = 2.0;

writeln("Domain (t = 0):");
writeln(u);

// apply the differential equation for nt iterations
var un : [0..#nx] real;

for n in 0..nt {
    u <=> un;

    foreach i in 1..<(nx-1) {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]);
    }
}

writeln("Domain (t = ", nt * dt,"):");
writeln(u);
