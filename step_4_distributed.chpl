use StencilDist;

// ----------------------------------------------------------------------------
// define a default setup
// ----------------------------------------------------------------------------
config const nx = 41;
const dx : real = 2.0 * pi / (nx - 1);
config const nt = 25;
config const nu = 0.3;
config const sigma = 0.2;
const dt = sigma * dx**2 / nu;


writeln("Running 1D Distributed Burger Simulation over: ");
writeln();
writeln("*--------(", nx, "x)--------* \t (dx = ", dx, ")");
writeln("0 \t\t  ", dx * nx);
writeln();
writeln("for ", nt * dt, " seconds (dt = ", dt, ")");
writeln("with: nu = ", nu);
writeln();

// ----------------------------------------------------------------------------
// create a distributed array to represent the computational Domain
// ----------------------------------------------------------------------------

// the one dimensional range of discrete points that we are interested in
const SpaceInner = {1..nx};
const SpaceFull = {0..nx+1};

// a distributed map which optimized for stencil operations
const CompDom = SpaceInner dmapped Stencil(
        SpaceInner,
        fluff=(1,),     // each local only needs to know about 1 adjacent point
        periodic=true   // we want to use periodic boundary conditions
    );

const FullDom = SpaceFull dmapped Stencil(
        SpaceInner,
        fluff=(1,),     // each local only needs to know about 1 adjacent point
        periodic=true   // we want to use periodic boundary conditions
    );

// the distributed array (defined over the full range)
var u : [FullDom] real;

// ----------------------------------------------------------------------------
// set up the initial conditions
// ----------------------------------------------------------------------------

u = 1.0;
u[(0.5 / dx):int..(1.0 / dx + 1):int] = 2.0;
u.updateFluff();

writeln("Domain (t = 0):");
writeln(u[CompDom]);

// ----------------------------------------------------------------------------
// apply the differential equation for nt iterations
// ----------------------------------------------------------------------------

var un : [FullDom] real;

for n in 0..nt {
    u <=> un;
    un.updateFluff();

    forall i in CompDom {
        u[i] = un[i] - un[i] * dt / dx * (un[i] - un[i - 1]) + nu * dt / dx**2 *
            (un[i + 1] - 2 * un[i] + un[i - 1]);
    }

    // cyclical boundary conditions
    u[0] = un[0] - un[0] * dt / dx * (un[0] - un[nx - 2]) + nu * dt / dx**2 *
        (un[1] -  2 * un[0] + un[nx - 2]);
    u[nx] = u[0];
}

writeln("Domain (t = ", nt * dt,"):");
writeln(u[CompDom]);
