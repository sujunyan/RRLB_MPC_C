#ifndef __PEMPC_H__
#define __PEMPC_H__

#include "MPC_problem.h"
#include <stddef.h>

#define USE_OPENMP (1)


#define y_size (N*(nx+nu)+nx)
#define lam_size ((N+1)*nx)

// self-defined linear algebra routine -----------------------
void pempc_aAxpy(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y);

void pempc_aAx(const char trans, const double a, const double* A, 
               const double* x, const size_t M, const size_t N, double* y);

void pempc_axpby(const double a, const double *x, const double b, 
                const double *y, const size_t N, double * out);

// main ALADIN algorithm -------------------------
/**
 * @brief Solve the coupled QP problem using Riccati backward-forward sweep.
 * @param xi1 Decoupled state sequence [x1, ..., xN].
 * @param xi2 Decoupled control sequence [u0, ..., u_{N-1}].
 * @param z1 Reference state sequence.
 * @param z2 Reference control sequence.
 * @param x0 Initial state vector of size nx.
 * @param z_out Output combined primal variable [x1...xN, u0...uN-1].
 * @param lam_inc_out Output dual variable increment [lam1...lamN].
 */
void solve_couple(const double *xi1, const double *xi2, const double *z1, const double *z2, 
                  const double *x0, double *z_out, double *lam_inc_out);
/**
 * @brief Solve the first ALADIN subproblem for all stages in parallel.
 * @param z1_in Input reference primal variable of size N*nx.
 * @param lam_in Input local dual variable of size N*nx.
 * @param xi1_out Output optimized primal variable of size N*nx.
 */
void solve_decouple1(const double* z1_in, const double* lam_in, double* xi1_out);
/**
 * @brief Solve the second ALADIN subproblem for all stages in parallel.
 * @param z2_in Input reference primal variable of size N*nu.
 * @param lam_in Input local dual variable of size N*nx.
 * @param xi2_out Output optimized primal variable of size N*nu.
 */
void solve_decouple2(const double* z2_in, const double* lam_in, double* xi2_out);

void mpc_get_control(const double *x0, size_t max_iter, double tol, 
            const double *zin ,const double *lam_in, 
            double *zout, double *lam_out, double *u0);

// MPT related functions --------------------------------
unsigned long mpt_eval(const double *X, const size_t i_start, double *U, 
const size_t domain, const size_t range, const size_t mpt_nr, const int * mpt_nc, 
const double *A, const double *B, const double *F, const double *G, 
const double *HTB, const double *GTB, const double *FTB);


#endif