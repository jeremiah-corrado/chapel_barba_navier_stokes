use util;

// define default simulation parameters
config const nx = 41,
             nt = 25,
             dt = 0.025;

const dx : real = 2.0 / (nx - 1);

writeln("Running 1D Non-Linear Convection Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t   ", dx * (nx - 1));
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");

// create an array to represent the computational Domain
var u : [0..<nx] real;

// set up the initial conditions
u = 1.0;
u[(0.5 / dx):int..<(1.0 / dx + 1):int] = 2.0;

writeln("\nDomain (t = 0):");
writeln(u);

// apply the fd equation for nt iterations
var un = u;
for n in 0..#nt {
    u <=> un;

    forall i in 1..<nx {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i-1]);
    }
}

writeln("\nDomain (t = ", nt * dt,"):");
writeln(u);

write_array_to_file("./sim_output/step_2/ch_u.txt", u);
write_array_to_file("./sim_output/step_2/ch_x.txt", linspace(0.0, 2.0, nx));
