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

const cdom = {0..<nx, 0..<ny};
const cdom_inner: subdomain(cdom) = cdom.expand((-1, -1));

var p : [cdom] real = 0.0; // pressure scalar
var u : [cdom] real = 0.0; // x component of momentum
var v : [cdom] real = 0.0; // y component of momentum

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

    uv_boundary(u, v);
    p_boundary(p);

    var b : [cdom] real;

    // run simulation for nt time steps
    for t_step in 0..#nt {
        // un <=> u;
        // vn <=> v;
        moveSwap(u, un);
        moveSwap(v, vn);

        // solve for the component of p that depends solely on u and v
        comp_b(b, un, vn);

        // iteratively solve for pressure
        for iteration in 0..#nit {
            // p <=> pn;
            moveSwap(p, pn);
            p_np1(p, pn, b);
            p_boundary(p);
        }

        // solve for u and v using the updated pressure values
        u_np1(u, un, vn, p);
        v_np1(v, un, vn, p);

        // apply boundary conditions to u and v
        uv_boundary(u, v);
    }
}

proc u_np1(ref u, const ref un, const ref vn, const ref p) {
    foreach (i, j) in cdom_inner {
        u[i, j] = un[i, j] -
            (un[i, j] * (dt / dx) * (un[i, j] - un[i, j-1])) -
            (vn[i, j] * (dt / dy) * (un[i, j] - un[i-1, j])) -
            (dt / (rho**2 * dx) * (p[i+1, j] - p[i-1, j])) +
            nu * (
                (dt / dx**2) * (un[i+1, j] - 2.0 * un[i, j] + un[i-1, j]) +
                (dt / dy**2) * (un[i, j+1] - 2.0 * un[i, j] + un[i, j-1])
            );
    }
}

proc v_np1(ref v, const ref un, const ref vn, const ref p) {
    foreach (i, j) in cdom_inner {
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

// compute p_n+1 from p_n and b
proc p_np1(ref p, const ref pn, const ref b) {
    foreach (i, j) in cdom_inner {
        p[i, j] = (
                    dy**2 * (pn[i+1, j] - pn[i-1, j]) +
                    dx**2 * (pn[i, j+1] - pn[i, j-1])
                ) / dxy2 - b[i, j];
    }
}

// compute b from u and v
proc comp_b(ref b, const ref u, const ref v) {
    var du, dv: real;
    foreach (i, j) in cdom_inner {
        du = u[i+1,j] - u[i-1,j];
        dv = v[i,j+1] - v[i,j-1];

        b[i, j] =
            (1.0 / dt) * (du / (2.0 * dx) + dv / (2.0 * dy)) -
            (du * du / (4.0 * dx**2)) -
            (dv * dv / (4.0 * dy**2)) -
            (
                (u[i,j+1] - u[i,j-1]) * (v[i+1,j] - v[i-1,j]) /
                (2.0 * dy * dx)
            );
    }

    b *= rho * dx**2 * dy**2 / dxy2;
}
proc uv_boundary(ref u, ref v) {
    u[.., 0] = 0.0;
    u[0, ..] = 0.0;
    u[.., ny - 1] = 1.0;
    u[nx - 1, ..] = 0.0;

    v[.., 0] = 0.0;
    v[0, ..] = 0.0;
    v[.., ny - 1] = 0.0;
    v[nx - 1, ..] = 0.0;
}

proc p_boundary(ref p) {
    p[.., 0] = p[.., 1];
    p[0, ..] = p[1, ..];
    p[.., ny - 1] = 0.0;
    p[nx - 1, ..] = p[nx - 2, ..];
}
