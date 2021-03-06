use util;
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

const cdom = {0..<nx, 0..<ny};
const cdom_inner: subdomain(cdom) = cdom.expand((-1, -1));

var p : [cdom] real = 1.0; // pressure scalar
var u : [cdom] real = 0.0; // x component of momentum
var v : [cdom] real = 0.0; // y component of momentum

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
    // while iteration <= 3 {
        moveSwap(u, un);
        moveSwap(v, vn);

        // compute the portion of p that depends only on u and v
        comp_b(b, u, v);

        // iteratively solve for p
        for p_iter in 0..#nit {
            moveSwap(p, pn);
            p_np1(p, pn, b);
            p_boundary(p);
        }

        // compute u and v
        u_np1(u, un, vn, p);
        v_np1(v, un, vn, p);
        u_boundary(u);
        v_boundary(v);

        // compute the relative change in u (have we reached steady state yet?)
        udiff = ((+ reduce u) - (+ reduce un)) / (+ reduce u);

        writeln("iteration: ", iteration, " udiff: ", udiff);
        iteration += 1;
    }
}

proc comp_b(ref b, const ref u, const ref v) {
    var du, dv: real;
    foreach ((i_m, i, i_p), (j_m, j, j_p)) in x_cyclical(cdom_inner) {
        du = u[i, j_p] - u[i, j_m];
        dv = v[i_p, j] - v[i_m, j];

        b[i, j] =
            (1.0 / dt) * (du / (2.0 * dx) + dv / (2.0 * dy)) -
            (du * du / (4.0 * dx**2)) -
            (dv * dv / (4.0 * dy**2)) -
            (
                (u[i_p, j] - u[i_m, j]) * (v[i, j_p] - v[i, j_m]) /
                (2.0 * dy * dx)
            );
    }

    b *= rho * dx**2 * dy**2 / dxy2;
}

proc p_np1(ref p, const ref pn, const ref b) {
    foreach ((i_m, i, i_p), (j_m, j, j_p)) in x_cyclical(cdom_inner) {
        p[i, j] = (
                    dy**2 * (pn[i, j_p] - pn[i, j_m]) +
                    dx**2 * (pn[i_p, j] - pn[i_m, j])
                ) / dxy2 - b[i, j];
    }
}

proc u_np1(ref u, const ref un, const ref vn, const ref p) {
    foreach ((i_m, i, i_p), (j_m, j, j_p)) in x_cyclical(cdom_inner) {
        u[i, j] = un[i, j] -
            (un[i, j] * (dt / dx) * (un[i, j] - un[i, j_m])) -
            (vn[i, j] * (dt / dy) * (un[i, j] - un[i_m, j])) -
            (dt / (rho**2 * dx) * (p[i_p, j] - p[i_m, j])) +
            nu * (
                (dt / dx**2) * (un[i_p, j] - 2.0 * un[i, j] + un[i_m, j]) +
                (dt / dy**2) * (un[i, j_p] - 2.0 * un[i, j] + un[i, j_m])
            ) +
            F * dt;
    }
}

proc v_np1(ref v, const ref un, const ref vn, const ref p) {
    foreach ((i_m, i, i_p), (j_m, j, j_p)) in x_cyclical(cdom_inner)  {
        v[i, j] = vn[i, j] -
            un[i, j] * (dt / dx) * (vn[i, j] - vn[i, j_m]) -
            vn[i, j] * (dt / dy) * (vn[i, j] - vn[i_m, j]) -
            dt / (rho**2 * dy) * (p[i, j_p] - p[i, j_m]) +
            nu * (
                (dt / dx**2) * (vn[i_p, j] - 2.0 * vn[i, j] + vn[i_m, j]) +
                (dt / dy**2) * (vn[i, j_p] - 2.0 * vn[i, j] + vn[i, j_m])
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

// a custom iterator to handle the cyclical boundary on the left and right walls
iter x_cyclical(d_inner) {
    // // left wall (x = 0)
    // for j in d_inner.dim(1) {
    //     yield ((nx - 1, 0, 1), j);
    // }
    // // inner domain
    // for (i, j) in d_inner {
    //     yield ((i - 1, i, i + 1), j);
    // }
    // // right wall (x = 2)
    // for j in d_inner.dim(1) {
    //     yield ((nx - 2, nx - 1, 0), j);
    // }


    // left wall (x = 0)
    for i in d_inner.dim(0) {
        yield ((i - 1, i, i + 1), (ny - 1, 0, 1));
    }
    // inner domain
    for (i, j) in d_inner {
        yield ((i - 1, i, i + 1), (j - 1, j, j + 1));
    }
    // right wall (x = 2)
    for i in d_inner.dim(0) {
        yield ((i - 1, i, i + 1), (ny - 2, ny - 1, 0));
    }
}
