use util;

// define u-function directly
var ufunc = lambda(t:real, x:real, nu:real) {
    return -2*nu*(-(-8*t + 2*x)*exp(-(-4*t + x)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)) -
        (-8*t + 2*x - 4*pi)*exp(-(-4*t + x - 2*pi)**2/(4*nu*(t + 1)))/(4*nu*(t + 1)))/(exp(-(-4*t + x - 2*pi)**2/(4*nu*(t + 1))) +
        exp(-(-4*t + x)**2/(4*nu*(t + 1)))) + 4;
};

// define default simulation parameters
config const nx = 101,
             nt = 100,
             nu = 0.07;

const dx = 2 * pi / (nx - 1),
      dt = dx * nu;

writeln("Running 1D Diffusion Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t  ", dx * (nx - 1));
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");

// setup initial conditions
const x = linspace(0.0, 2 * pi, nx);
const cdom = 0..<nx;
const cdom_inner = 1..<(nx-1);

var u : [cdom] real;
[i in cdom] u[i] = ufunc(0, x[i], nu);

writeln("u(t = 0, x):");
writeln(u);

// apply the fd equation for nt iterations
var un = u;
for i in 0..#nt {
    u <=> un;

    forall i in cdom_inner {
        u[i] = un[i] - un[i] * dt / dx *(un[i] - un[i-1]) + nu * dt / dx**2 *
                (un[i+1] - 2 * un[i] + un[i-1]);
    }

    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx-2]) + nu * dt / dx**2 *
                (un[1] - 2 * un[0] + un[nx-2]);
    u[nx-1] = u[0];
}

writeln("\nDomain (t = ", nt * dt,", x):");
writeln(u);

write_array_to_file("./sim_output/step_4/ch_u.txt", u);
write_array_to_file("./sim_output/step_4/ch_x.txt", x);
