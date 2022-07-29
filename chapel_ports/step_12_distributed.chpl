use util;
import Memory.Initialization.moveSwap;

config const nt = 100; // number of time steps
config const dt = 0.001; // temporal resolution
config const nit = 100; // number of diffusion resolution iterations

config const nx = 41; // x spatial-resolution
config const ny = 41; // y spatial-resolution
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
const dxy2 = 2.0 * (dx**2 + dy**2);

config const rho = 1;
config const nu = 0.1;
config const F = 1;

const cdom = {0..<nx, 0..<ny};
const cdom_inner: subdomain(cdom) = cdom.expand((-1, -1));
const CDOM = cdom dmapped Stencil(cdom_inner, fluff=(1,1), periodic=true);

var p : [cdom] real; // pressure scalar
var u : [cdom] real; // x component of momentum
var v : [cdom] real; // y component of momentum

channel_flow_sim(u, v, p, 0.001);

write_array_to_file("./sim_output/step_12/ch_u.txt", u);
write_array_to_file("./sim_output/step_12/ch_v.txt", v);
write_array_to_file("./sim_output/step_12/ch_p.txt", p);
write_array_to_file("./sim_output/step_12/ch_x.txt", linspace(0.0, 2.0, nx));
write_array_to_file("./sim_output/step_12/ch_y.txt", linspace(0.0, 2.0, ny));

proc channel_flow_sim(ref u, ref v, ref p, udiff_thresh: real) {
    var udiff = 1.0;
    var iteration = 0;

    var un = u;
    var vn = v;
    var pn = p;

    var b : [cdom] real;

    while udiff > udiff_thresh {
        moveSwap(u, un);
        moveSwap(v, vn);

        comp_b(b, u, v);

        // iteratively solve for pressure
        for iteration in 0..#nit {
            moveSwap(p, pn);
            p_np1(p, pn, b);
            p_boundary(p);
        }

        u_np1(u, un, vn, p);
        v_np1(v, un, vn, p);

        u_boundary(u);
        v_boundary(v);

        udiff = ((+ reduce u) - (+ reduce un)) / (+ reduce u);
        iteration += 1;

        writeln("iteration: ", iteration, " udiff: ", udiff);
    }
}

proc comp_b(ref b, const ref u, const ref v) {
    var du, dv: real;
    foreach ((i_m, i, i_p), j) in x_cyclical(cdom_inner) {
        du = u[i_p, j] - u[i_m, j];
        dv = v[i, j+1] - v[i, j-1];

        b[i, j] =
            (1.0 / dt) * (du / (2.0 * dx) + dv / (2.0 * dy)) -
            (du * du / (4.0 * dx**2)) -
            (dv * dv / (4.0 * dy**2)) -
            (
                (u[i, j+1] - u[i, j-1]) * (v[i_p, j] - v[i_m, j]) /
                (2.0 * dy * dx)
            );
    }

    b *= rho * dx**2 * dy**2 / dxy2;
}

proc p_np1(ref p, const ref pn, const ref b) {
    foreach ((i_m, i, i_p), j) in x_cyclical(cdom_inner) {
        p[i, j] = (
                    dy**2 * (pn[i_p, j] - pn[i_m, j]) +
                    dx**2 * (pn[i, j+1] - pn[i, j-1])
                ) / dxy2 - b[i, j];
    }
}

proc u_np1(ref u, const ref un, const ref vn, const ref p) {
    foreach ((i_m, i, i_p), j) in x_cyclical(cdom_inner) {
        u[i, j] = un[i, j] -
            (un[i, j] * (dt / dx) * (un[i, j] - un[i, j-1])) -
            (vn[i, j] * (dt / dy) * (un[i, j] - un[i_m, j])) -
            (dt / (rho**2 * dx) * (p[i_p, j] - p[i_m, j])) +
            nu * (
                (dt / dx**2) * (un[i_p, j] - 2.0 * un[i, j] + un[i_m, j]) +
                (dt / dy**2) * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])
            ) +
            F * dt;
    }
}

proc v_np1(ref v, const ref un, const ref vn, const ref p) {
    foreach ((i_m, i, i_p), j) in x_cyclical(cdom_inner)  {
        v[i, j] = vn[i, j] -
            un[i, j] * (dt / dx) * (vn[i, j] - vn[i, j-1]) -
            vn[i, j] * (dt / dy) * (vn[i, j] - vn[i_m, j]) -
            dt / (rho**2 * dy) * (p[i, j+1] - p[i, j-1]) +
            nu * (
                (dt / dx**2) * (vn[i_p, j] - 2.0 * vn[i, j] + vn[i_m, j]) +
                (dt / dy**2) * (vn[i, j+1] - 2.0 * vn[i, j] + vn[i, j-1])
            );
    }
}

proc u_boundary(ref u) {
    u[.., 0] = 0.0;
    u[.., ny-1] = 0.0;
}

proc v_boundary(ref v) {
    v[.., 0] = 0.0;
    v[.., ny-1] = 0.0;
}

proc p_boundary(ref p) {
    p[.., 0] = p[.., 1];
    p[.., ny-1] = p[.., ny-2];
}

// a custom iterator to handle the cyclical boundary on the left and right walls
iter x_cyclical(param tag: iterKind, A) where tag == iterKind.standalone {
    coforall loc in A.targetLocales on loc {
        for (i, j) in  A.localSubdomain(here) {
            yield ((i - 1, i, i + 1), j);
        }
    }
    // left wall (x = 0)
    for j in d_inner.dim(1) {
        yield ((nx -1, 0, 1), j);
    }
    // right wall (x = 2)
    for j in d_inner.dim(1) {
        yield ((nx - 2, nx - 1, 0), j);
    }
}
