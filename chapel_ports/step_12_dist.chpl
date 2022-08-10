use util;
use StencilDist;
import Memory.Initialization.moveSwap;

config const nt = 10; // number of time steps
config const dt = 0.01; // temporal resolution
config const nit = 50; // number of diffusion resolution iterations

config const nx = 41; // x spatial-resolution
config const ny = 41; // y spatial-resolution
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
const dxy2 = 2.0 * (dx**2 + dy**2);

config const rho = 1;
config const nu = 0.1;
config const F = 1;

config const write_data = false;

// Define the usual domain with an extra 2-element buffer in the x-domain
//  this is done so that the periodicity in the x-direction can be ignored.
const cdom = {0..<(nx+2), 0..<ny};
const cdom_inner: subdomain(cdom) = cdom.expand((-2, -1)); // the region over which we'll do stencil computations
const cdom_actual: subdomain(cdom) = cdom.expand((-1, 0)); // the region that contains sensible solution information

// define the distributed domain with cdom_inner as the bounding box, and periodicity active
const CDOM = cdom dmapped Stencil(cdom_inner, fluff=(1,1), periodic=true);
const CDOM_INNNER = CDOM[cdom_inner];

var p : [CDOM] real = 0.0; // pressure scalar
var u : [CDOM] real = 0.0; // x component of momentum
var v : [CDOM] real = 0.0; // y component of momentum

channel_flow_sim(u, v, p, 0.001);

if write_data {
    write_array_to_file("./sim_output/step_12/ch_u.txt", u[cdom_actual]);
    write_array_to_file("./sim_output/step_12/ch_v.txt", v[cdom_actual]);
    write_array_to_file("./sim_output/step_12/ch_p.txt", p[cdom_actual]);
    write_array_to_file("./sim_output/step_12/ch_x.txt", linspace(0.0, 2.0, nx));
    write_array_to_file("./sim_output/step_12/ch_y.txt", linspace(0.0, 2.0, ny));
}

proc channel_flow_sim(ref u, ref v, ref p, udiff_thresh: real) {
    var udiff = 1.0;
    var i = 0;

    var un : [CDOM] real = u;
    var vn : [CDOM] real = v;
    var pn : [CDOM] real = p;

    var b : [CDOM] real;

    while udiff > udiff_thresh {
    // while iteration <= 3 {
        moveSwap(u, un);
        moveSwap(v, vn);

        // compute the portion of p that depends only on u and v
        comp_b(b, u, v);

        b.updateFluff();

        // iteratively solve for p
        for p_iter in 0..#nit {
            moveSwap(p, pn);
            p_np1(p, pn, b);
            p_boundary(p);
            p.updateFluff();
        }

        // compute u and v
        cobegin {
            {
                u_np1(u, un, vn, p);
                u_boundary(u);
                u.updateFluff();
            }
            {
                v_np1(v, un, vn, p);
                v_boundary(v);
                v.updateFluff();
            }
        }

        // compute the relative change in u (have we reached steady state yet?)
        udiff = ((+ reduce u) - (+ reduce un)) / (+ reduce u);
        i += 1;

        // writeln("iteration: ", i, " udiff: ", udiff);
    }

    writeln("ran for ", i, " iterations (final udiff = ", udiff, ")");
}

proc comp_b(ref b : [] real, const ref u, const ref v) {
    forall (i, j) in CDOM_INNNER with (var du: real, var dv: real) {
        du = u[i, j+1] - u[i, j-1];
        dv = v[i+1, j] - v[i-1, j];

        b[i, j] = rho * (1.0 / dt) *
            (du / (2.0 * dx) + dv / (2.0 * dy)) -
            (du / (2.0 * dx))**2 -
            (dv / (2.0 * dy))**2 -
            2.0 * (
                (u[i+1, j] - u[i-1, j]) / (2.0 * dy) *
                (v[i, j+1] - v[i, j-1]) / (2.0 * dx)
            );
    }
}

proc p_np1(ref p : [] real, const ref pn, const ref b) {
    forall (i, j) in CDOM_INNNER {
        p[i, j] = (
                    dy**2 * (pn[i, j+1] + pn[i, j-1]) +
                    dx**2 * (pn[i+1, j] + pn[i-1, j])
                ) / dxy2 - dx**2 * dy**2 / dxy2 * b[i, j];
    }
}

proc u_np1(ref u : [] real, const ref un, const ref vn, const ref p) {
    forall (i, j) in CDOM_INNNER {
        u[i, j] = un[i, j] -
            (un[i, j] * (dt / dx) * (un[i, j] - un[i, j-1])) -
            (vn[i, j] * (dt / dy) * (un[i, j] - un[i-1, j])) -
            (dt / (rho**2 * dx) * (p[i+1, j] - p[i-1, j])) +
            nu * (
                (dt / dx**2) * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) +
                (dt / dy**2) * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])
            ) +
            F * dt;
    }
}

proc v_np1(ref v : [] real, const ref un, const ref vn, const ref p) {
    forall (i, j) in CDOM_INNNER  {
        v[i, j] = vn[i, j] -
            un[i, j] * (dt / dx) * (vn[i, j] - vn[i, j-1]) -
            vn[i, j] * (dt / dy) * (vn[i, j] - vn[i-1, j]) -
            dt / (rho**2 * dy) * (p[i, j+1] - p[i, j-1]) +
            nu * (
                (dt / dx**2) * (vn[i+1, j] - 2.0 * vn[i, j] + vn[i-1, j]) +
                (dt / dy**2) * (vn[i, j+1] - 2.0 * vn[i, j] + vn[i, j-1])
            );
    }
}

proc u_boundary(ref u) {
    u[0, ..] = 0.0;
    u[nx-1, ..] = 0.0;
}

proc v_boundary(ref v) {
    v[0, ..] = 0.0;
    v[nx-1, ..] = 0.0;
}

proc p_boundary(ref p) {
    p[0, ..] = p[1, ..];
    p[nx-1, ..] = p[nx-2, ..];
}
