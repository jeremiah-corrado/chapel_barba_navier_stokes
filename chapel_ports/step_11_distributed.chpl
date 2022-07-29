use stencilDist;
use util;

config const nt = 100; // number of time steps
config const dt = 0.001; // temporal resolution
config const nit = 100; // number of diffusion resolution iterations

config const nx = 41; // x spatial-resolution
config const ny = 41; // y spatial-resolution
const dx = 2.0 / (nx - 1);
const dy = 2.0 / (ny - 1);
const dxy2 = 2.0 * (dx**2 + dy**2);

config const rho = 1;
config const nu = 1;

const cdom = {0..<nx, 0..<ny};
const cdom_inner: subdomain(cdom) = cdom.expand((-1, -1));
const CDOM = cdom dmapped Stencil(cdom_inner, fluff=(1,1));

var p : [CDOM] real; // pressure scalar
var u : [CDOM] real; // x component of flow
var v : [CDOM] real; // y component of flow

cavity_flow_sim(u, v, p);

write_array_to_file("./sim_output/step_11/ch_u.txt", u);
write_array_to_file("./sim_output/step_11/ch_v.txt", v);
write_array_to_file("./sim_output/step_11/ch_p.txt", p);
write_array_to_file("./sim_output/step_11/ch_x.txt", linspace(0.0, 2.0, nx));
write_array_to_file("./sim_output/step_11/ch_y.txt", linspace(0.0, 2.0, ny));

proc cavity_flow_sim(ref u, ref v, ref p) {
    // temporary copies of computational domain
    var un = u;
    var vn = v;
    var pn = p;

    var b : [cdom] real;

    // run simulation for nt time steps
    for t_step in 0..#nt {
        u <=> un;
        v <=> vn;

        // solve for the component of p that depends solely on u and v
        comp_b(b, un, vn);

        // iteratively solve for pressure
        for iteration in 0..#nit {
            p <=> pn;
            p_np1(p, pn, b);
        }

        // solve for u and v using the updated pressure values
        cobegin {
            u_np1(u, un, vn, p);
            v_np1(v, un, vn, p);
        }

        // apply boundary conditions to u and v
        apply_boundary(u, v);

        u.updateFluff();
        v.updateFluff();
    }
}

proc u_np1(ref u, const ref un, const ref vn, const ref p) {
    forall (i, j) in cdom_inner {
        u[i, j] = un[i, j] -
            un[i, j] * (dt / dx) * (un[i, j] - un[i-1, j]) -
            vn[i, j] * (dt / dy) * (un[i, j] - un[i, j-1]) -
            dt / (rho**2 * dx) * (p[i+1, j] - p[i-1, j]) +
            nu * (
                (dt / dx**2) * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) +
                (dt / dy**2) * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])
            );
    }
}

proc v_np1(ref v, const ref un, const ref vn, const ref p) {
    forall (i, j) in cdom_inner {
        v[i, j] = vn[i, j] -
            un[i, j] * (dt / dx) * (vn[i, j] - vn[i-1, j]) -
            vn[i, j] * (dt / dy) * (vn[i, j] - vn[i, j-1]) -
            dt / (rho**2 * dx) * (p[i, j+1] - p[i, j-1]) +
            nu * (
                (dt / dx**2) * (vn[i+1, j] - 2.0 * vn[i, j] + vn[i-1, j]) +
                (dt / dy**2) * (vn[i, j+1] - 2.0 * vn[i, j] + vn[i, j-1])
            );
    }
}

// compute p_n+1 from p_n and b
proc p_np1(ref p, const ref pn, const ref b) {
    forall (i, j) in cdom_inner {
        p[i, j] = (
                    dy**2 * (pn[i+1, j] - pn[i-1, j]) +
                    dx**2 * (pn[i, j+1] - pn[i, j-1])
                ) / dxy2 * b[i, j];
    }
}

// compute b from u and v
proc comp_b(ref b, const ref u, const ref v) {
    forall (i, j) in cdom_inner with (var du: real, var dv: real) {
        du = u[i+1,j] - u[i-1,j];
        dv = v[i,j+1] - v[i,j-1];

        b[i, j] = ((1.0 / dt) * du / (2.0 * dx) + dv / (2.0 * dy)) -
            ((0.25 / dx**2) * du * du) -
            ((0.25 / dy**2) * dv * dv) -
            (0.5 / (dx * dy)) * (u[i,j+1] - u[i,j-1]) * (v[i+1,j] - v[i-1,j]);
    }

    b *= rho * dx**2 * dy**2 / dxy2;
}

proc apply_boundary(ref u, ref v) {
    u[.., 0] = 0.0;
    u[0, ..] = 0.0;
    u[.., ny - 1] = 1.0;
    u[nx - 1, ..] = 0.0;

    v[.., 0] = 0.0;
    v[0, ..] = 0.0;
    v[.., ny - 1] = 0.0;
    v[nx - 1, ..] = 0.0;
}
