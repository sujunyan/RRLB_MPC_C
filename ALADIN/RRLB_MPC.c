#include <string.h>
#include <stdio.h>
#include <math.h>
#include "RRLB_MPC.h"
#include <time.h>
#include <float.h>
#include "mpQP_data.h"
#include <assert.h>
#if USE_OPENMP
#include <omp.h>
#endif

const double MPT_ABSTOL = 1.000000e-8;

/**
 * @brief mpt evaluation function using sequential search
 * @param  *X :   nux1 vector, the parameter theta_k
 * @param  *U :   nux1 vector, the output optimal solution of the mpqp
 * @param  istart: the start index for searching region, to speed up the search. Normally set to 0.
*/
static unsigned long MPT_func_k_( double *X, double *U, const size_t istart){
    return mpt_eval(X,istart,U,nu,nu,kMPT_NR,kMPT_NC,kMPT_A,kMPT_B,kMPT_F,kMPT_G,kMPT_HTB,kMPT_GTB,kMPT_FTB);
}

/**
 * @brief Solve a linear system Ax = b for n x n matrix using Gaussian elimination with partial pivoting.
 */
static void solve_linear_system(const double *A, const double *b, size_t n, double *x) {
    double mat[n][n + 1];
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            mat[i][j] = A[i + j * n];
        }
        mat[i][n] = b[i];
    }
    for (size_t i = 0; i < n; i++) {
        size_t max = i;
        for (size_t j = i + 1; j < n; j++) {
            if (fabs(mat[j][i]) > fabs(mat[max][i])) max = j;
        }
        for (size_t j = i; j <= n; j++) {
            double tmp = mat[i][j];
            mat[i][j] = mat[max][j];
            mat[max][j] = tmp;
        }
        for (size_t j = i + 1; j < n; j++) {
            double factor = mat[j][i] / mat[i][i];
            for (size_t k = i; k <= n; k++) {
                mat[j][k] -= factor * mat[i][k];
            }
        }
    }
    for (int i = (int)n - 1; i >= 0; i--) {
        double sum = 0;
        for (size_t j = i + 1; j < n; j++) {
            sum += mat[i][j] * x[j];
        }
        x[i] = (mat[i][n] - sum) / mat[i][i];
    }
}

static double relax_barrier_d1_gen(double delta, double x) {
    if (x >= delta) return -1.0 / x;
    else return (x - 2.0 * delta) / (delta * delta);
}

static double relax_barrier_d2_gen(double delta, double x) {
    if (x >= delta) return 1.0 / (x * x);
    else return 1.0 / (delta * delta);
}

static void f1_ALADIN_gradient_single_c(const double *xi1_j, const double *z1_j, const double *lam_j, int terminal_flag, double *g_j) {
    const double *Q_mat = terminal_flag ? P : Q;
    pempc_aAx('n', 2.0, Q_mat, xi1_j, nx, nx, g_j);
    for (size_t i = 0; i < nx; i++) g_j[i] += q[i];
    if (!terminal_flag) {
        for (size_t r = 0; r < mx; r++) {
            double Cx_r_xi = 0;
            for (size_t i = 0; i < nx; i++) Cx_r_xi += Cx[r + i * mx] * xi1_j[i];
            double val = dx[r] - Cx_r_xi;
            double rlb_d1 = relax_barrier_d1_gen(delta_relax, val);
            double factor = rho * wx[r] * rlb_d1;
            for (size_t i = 0; i < nx; i++) g_j[i] -= factor * Cx[r + i * mx];
        }
    }
    for (size_t i = 0; i < nx; i++) g_j[i] += lam_j[i];
    double xi_diff[nx];
    for (size_t i = 0; i < nx; i++) xi_diff[i] = xi1_j[i] - z1_j[i];
    pempc_aAxpy('n', 1.0, ALADIN_Q_bar, xi_diff, nx, nx, g_j);
}

static void f1_ALADIN_Hessian_single_c(const double *xi1_j, int terminal_flag, double *H_j) {
    const double *Q_mat = terminal_flag ? P : Q;
    for (size_t i = 0; i < nx * nx; i++) H_j[i] = 2.0 * Q_mat[i];
    if (!terminal_flag) {
        for (size_t r = 0; r < mx; r++) {
            double Cx_r_xi = 0;
            for (size_t i = 0; i < nx; i++) Cx_r_xi += Cx[r + i * mx] * xi1_j[i];
            double val = dx[r] - Cx_r_xi;
            double rlb_d2 = relax_barrier_d2_gen(delta_relax, val);
            double factor = rho * wx[r] * rlb_d2;
            for (size_t i = 0; i < nx; i++) {
                for (size_t k = 0; k < nx; k++) {
                    H_j[i + k * nx] += factor * Cx[r + i * mx] * Cx[r + k * mx];
                }
            }
        }
    }
    for (size_t i = 0; i < nx * nx; i++) H_j[i] += ALADIN_Q_bar[i];
}

void solve_decouple1(const double* z1_in, const double* lam_in, double* xi1_out){
    const double grad_tol = 1e-6;
    const double alpha = 1.0;
    memcpy(xi1_out, z1_in, sizeof(double) * N * nx);
    for (size_t k = 0; k < 10000; k++) {
        double delta_xi1[N * nx];
        double max_delta = 0;
#if USE_OPENMP
#pragma omp parallel for reduction(max:max_delta)
#endif
        for (size_t j = 0; j < N; j++) {
            int terminal_flag = (j == N - 1);
            const double *xi1_j = &xi1_out[j * nx];
            const double *z1_j = &z1_in[j * nx];
            const double *lam_j = &lam_in[j * nx];
            double g_j[nx];
            f1_ALADIN_gradient_single_c(xi1_j, z1_j, lam_j, terminal_flag, g_j);
            double g_max = 0;
            for (size_t i = 0; i < nx; i++) {
                if (fabs(g_j[i]) > g_max) g_max = fabs(g_j[i]);
            }
            double delta_xi1_j[nx];
            if (g_max < grad_tol) {
                for (size_t i = 0; i < nx; i++) delta_xi1_j[i] = 0;
            } else {
                double H_j[nx * nx];
                f1_ALADIN_Hessian_single_c(xi1_j, terminal_flag, H_j);
                double minus_g_j[nx];
                for (size_t i = 0; i < nx; i++) minus_g_j[i] = -g_j[i];
                solve_linear_system(H_j, minus_g_j, nx, delta_xi1_j);
            }
            for (size_t i = 0; i < nx; i++) {
                delta_xi1[j * nx + i] = delta_xi1_j[i];
                double abs_d = fabs(delta_xi1_j[i]);
                if (abs_d > max_delta) max_delta = abs_d;
            }
        }
        if (max_delta < grad_tol) break;
        for (size_t i = 0; i < N * nx; i++) xi1_out[i] += alpha * delta_xi1[i];
    }
}

void solve_decouple2(const double* z2_in, const double* lam_in, double* xi2_out){
    double g[N * nu];
    double p_next[nx];
    memcpy(p_next, lam_in + (N-1)*nx, nx * sizeof(double));
    for (int k = (int)N - 1; k >= 0; k--) {
        pempc_aAx('t', 1.0, B, p_next, nx, nu, g + k*nu);
        if (k > 0) {
            double p_temp[nx];
            memcpy(p_temp, lam_in + (k-1)*nx, nx * sizeof(double));
            pempc_aAxpy('t', 1.0, A, p_next, nx, nx, p_temp);
            memcpy(p_next, p_temp, nx * sizeof(double));
        }
    }
#if USE_OPENMP
#pragma omp parallel for
#endif
    for (size_t k = 0; k < N; k++) {
        double theta_k[nu];
        const double *g_k = g + k*nu;
        const double *u_k = z2_in + k*nu;
        for (size_t i = 0; i < nu; i++) theta_k[i] = -g_k[i];
        pempc_aAxpy('n', -1.0, R, u_k, nu, nu, theta_k);
        MPT_func_k_(theta_k, xi2_out + k*nu, 0);
    }
}

void solve_couple(const double *xi1, const double *xi2, const double *z1, const double *z2, 
                  const double *x0, double *z_out, double *lam_inc_out){
    double p[(N+1)*nx];
    double l[N*nu];
    double qn[nx], sn[nu], u_tmp[nu], x_tmp[nx], xn[nx], un[nu];
    double ximz_tmp[nx];
    const double *xi1_N = xi1 + (N-1)*nx;
    const double *z1_N = z1 + (N-1)*nx;
    for (size_t i=0; i<nx; i++) ximz_tmp[i] = 2.0*xi1_N[i] - z1_N[i];
    pempc_aAx('n', -1.0, ALADIN_Q_bar, ximz_tmp, nx, nx, p + N*nx);
    for (int n = (int)N; n >= 1; n--) {
        const double *xi1_n_prev, *z1_n_prev;
        if (n == 1) {
            double zero_nx[nx];
	    memset(zero_nx, 0, sizeof(zero_nx));
            xi1_n_prev = zero_nx;
            z1_n_prev = zero_nx;
        } else {
            xi1_n_prev = xi1 + (n-2)*nx;
            z1_n_prev = z1 + (n-2)*nx;
        }
        const double *xi2_n = xi2 + (n-1)*nu;
        const double *z2_n = z2 + (n-1)*nu;
        for (size_t i=0; i<nx; i++) ximz_tmp[i] = 2.0*xi1_n_prev[i] - z1_n_prev[i];
        pempc_aAx('n', -1.0, ALADIN_Q_bar, ximz_tmp, nx, nx, qn);
        double uimz_n[nu];
        for (size_t i=0; i<nu; i++) uimz_n[i] = 2.0*xi2_n[i] - z2_n[i];
        pempc_aAx('n', -1.0, R, uimz_n, nu, nu, sn);
        const double *Ric_Lam_inv_n = Ric_Lam_inv + (n-1)*nu*nu;
        const double *Ric_L_n = Ric_L + (n-1)*nu*nx;
        double *ln = l + (n-1)*nu;
        double *pn = p + (n-1)*nx;
        double *p_next = p + n*nx;
        memcpy(u_tmp, sn, nu * sizeof(double));
        pempc_aAxpy('t', 1.0, B, p_next, nx, nu, u_tmp); 
        pempc_aAx('n', 1.0, Ric_Lam_inv_n, u_tmp, nu, nu, ln);
        pempc_aAxpy('t', 1.0, A, p_next, nx, nx, qn);
        pempc_aAxpy('t', -1.0, Ric_L_n, ln, nu, nx, qn);
        memcpy(pn, qn, nx * sizeof(double));
    }
    memcpy(xn, x0, nx * sizeof(double));
    for (size_t n = 1; n <= N; n++) {
        const double *Ric_Lam_inv_n = Ric_Lam_inv + (n-1)*nu*nu;
        const double *Ric_L_n = Ric_L + (n-1)*nu*nx;
        double *ln = l + (n-1)*nu;
        memcpy(u_tmp, ln, nu * sizeof(double));
        pempc_aAxpy('n', 1.0, Ric_L_n, xn, nu, nx, u_tmp);
        pempc_aAx('t', -1.0, Ric_Lam_inv_n, u_tmp, nu, nu, un);
        memcpy(z_out + N*nx + (n-1)*nu, un, nu * sizeof(double));
        pempc_aAx('n', 1.0, A, xn, nx, nx, x_tmp);
        pempc_aAxpy('n', 1.0, B, un, nx, nu, x_tmp);
        memcpy(xn, x_tmp, nx * sizeof(double));
        memcpy(z_out + (n-1)*nx, xn, nx * sizeof(double));
        const double *xi1_n = xi1 + (n-1)*nx;
        const double *z1_n = z1 + (n-1)*nx;
        for (size_t i=0; i<nx; i++) ximz_tmp[i] = xn[i] - (2.0*xi1_n[i] - z1_n[i]);
        pempc_aAx('n', -1.0, ALADIN_Q_bar, ximz_tmp, nx, nx, lam_inc_out + (n-1)*nx);
    }
}

void mpc_get_control(const double *x0_val, size_t max_iter, double tol, 
            const double *zin ,const double *lam_in, 
            double *zout, double *lam_out, double *u0){
    const size_t nx_total = N * nx;
    const size_t nu_total = N * nu;
    double xi1[N * nx];
    double xi2[N * nu];
    double z1[N * nx];
    double z2[N * nu];
    double lam[N * nx];
    double del_lam[N * nx];
    double z_combined[N * (nx + nu)];
    memcpy(z1, zin, nx_total * sizeof(double));
    memcpy(z2, zin + nx_total, nu_total * sizeof(double));
    memcpy(lam, lam_in, nx_total * sizeof(double));
    
    size_t m;
    for (m = 0; m < max_iter; m++){
        solve_decouple1(z1, lam, xi1);
        solve_decouple2(z2, lam, xi2);
        solve_couple(xi1, xi2, z1, z2, x0_val, z_combined, del_lam);
        
        for (size_t i = 0; i < nx_total; i++) lam[i] += del_lam[i];
        memcpy(z1, z_combined, nx_total * sizeof(double));
        memcpy(z2, z_combined + nx_total, nu_total * sizeof(double));

        double err = 0;
        for (size_t i = 0; i < nx_total; i++) {
            double d = fabs(z1[i] - xi1[i]);
            if (d > err) err = d;
        }
        for (size_t i = 0; i < nu_total; i++) {
            double d = fabs(z2[i] - xi2[i]);
            if (d > err) err = d;
        }
        if (err < tol) {
            break;
        }
    }
    // printf("Step: ALADIN converged in %zu iterations\n", m + 1);
    
    memcpy(u0, xi2, nu * sizeof(double));
    
    if (N > 1) {
        memcpy(zout, z1 + nx, (N - 1) * nx * sizeof(double));
    }
    memset(zout + (N - 1) * nx, 0, nx * sizeof(double));
    if (N > 1) {
        memcpy(zout + nx_total, z2 + nu, (N - 1) * nu * sizeof(double));
    }
    memset(zout + nx_total + (N - 1) * nu, 0, nu * sizeof(double));
    if (N > 1) {
        memcpy(lam_out, lam + nx, (N - 1) * nx * sizeof(double));
    }
    memset(lam_out + (N - 1) * nx, 0, nx * sizeof(double));
}

void pempc_axpby(const double a, const double *x, const double b, 
                const double *y, const size_t N_len, double * out){
    for (size_t i = 0; i<N_len; i++){
        out[i] = a*x[i] + b*y[i];
    }
}

void pempc_aAxpy(const char trans, const double a, const double* A_mat, 
               const double* x, const size_t M, const size_t N_len, double* y){
    for (size_t i = 0; i < (trans == 'n' ? M : N_len); i++) {
        double sum = 0;
        for (size_t j = 0; j < (trans == 'n' ? N_len : M); j++) {
            if (trans == 'n') sum += A_mat[i + j * M] * x[j];
            else sum += A_mat[j + i * M] * x[j];
        }
        y[i] += a * sum;
    }
}

void pempc_aAx(const char trans, const double a, const double* A_mat, 
               const double* x, const size_t M, const size_t N_len, double* y){
    size_t out_len = (trans == 'n' ? M : N_len);
    memset(y, 0, out_len * sizeof(double));
    pempc_aAxpy(trans, a, A_mat, x, M, N_len, y);
}

unsigned long mpt_eval(const double *X, const size_t i_start, double *U, 
const size_t domain, const size_t range, const size_t mpt_nr, const int * mpt_nc, 
const double *A_mpt, const double *B_mpt, const double *F_mpt, const double *G_mpt, 
const double *HTB, const double *GTB, const double *FTB){
    size_t abspos = 0;
    for (size_t ireg = 0; ireg < i_start; ireg++) abspos += mpt_nc[ireg];
    size_t iregmin = 0;
    unsigned long region = 0;
    for (size_t ireg0=0; ireg0<mpt_nr; ireg0++) {
        size_t ireg = (i_start + ireg0)%mpt_nr;
        int isinside = 1;
        int nc_curr = mpt_nc[ireg];
        for (int ic=0; ic<nc_curr; ic++) {
            double hx = 0;
            for (size_t ix=0; ix<domain; ix++) hx += A_mpt[abspos*domain+ic*domain+ix]*X[ix];
            if ((hx - B_mpt[abspos+ic]) > MPT_ABSTOL) {
                isinside = 0;
                break;
            } 
        }
        if (isinside == 1){
            region = ireg + 1;
            iregmin = ireg;
            break;
        }
        if (ireg == mpt_nr -1) abspos = 0;
        else abspos += mpt_nc[ireg];
    }
    for (size_t ix=0; ix<range; ix++) {
        double sx = 0;
        for (size_t jx=0; jx<domain; jx++) sx += F_mpt[iregmin*domain*range + ix*domain + jx]*X[jx];
        U[ix] = sx + G_mpt[iregmin*range + ix];
    }
    return region;
}
