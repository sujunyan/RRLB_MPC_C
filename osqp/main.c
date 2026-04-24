#include <stdio.h>
#include <string.h>
#include <time.h>
#include "workspace.h"
#include "osqp.h"

#include "MPC_problem.h"
// #include "qp_bound_data.h"
#include <stdlib.h> 

#define PRINT_DEBUG 1

/**
 * @brief Phase 3: Closed-loop simulation runner for OSQP.
 * This program simulates the plant dynamics using the OSQP-implemented MPC controller.
 * It outputs the state and control trajectories in CSV format for comparison with MATLAB.
 */

// Self-defined linear algebra routines (mimicking implementation in ALADIN)
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

void evaluate_cost_J(double* x, double* u, double* J) {
    double cost_x = 0.0;
    double Qx[nx];
    for (size_t i = 0; i < nx; i++) {
        Qx[i] = 0.0;
        for (size_t j = 0; j < nx; j++) {
            Qx[i] += Q[i * nx + j] * x[j];
        }
    }
    for (size_t i = 0; i < nx; i++) {
        cost_x += x[i] * Qx[i];
    }

    double cost_u = 0.0;
    double Ru[nu];
    for (size_t i = 0; i < nu; i++) {
        Ru[i] = 0.0;
        for (size_t j = 0; j < nu; j++) {
            Ru[i] += R[i * nu + j] * u[j];
        }
    }
    for (size_t i = 0; i < nu; i++) {
        cost_u += u[i] * Ru[i];
    }

    *J = cost_x + cost_u;
}


void update_x0(OSQPWorkspace *work, double* _x0){

    int m = work->data->m;

    double l_new[m];
    double u_new[m];

    double beq[nx*N];

    // Note that we cannot directly use work->data->l since OSQP have scaling internally.
    memcpy(l_new, QP_l, m * sizeof(double));
    memcpy(u_new, QP_u, m * sizeof(double));

    double current_xk[nx];
    memcpy(current_xk, _x0, nx * sizeof(double));

    for (size_t k = 0; k < N; k++) {
        double next_xk[nx];
        pempc_aAx('n', 1.0, A, current_xk, nx, nx, next_xk);
        memcpy(current_xk, next_xk, nx * sizeof(double));

        for (size_t i = 0; i < nx; i++) {
            beq[k * nx + i] = current_xk[i];
        }
    }
    memcpy(l_new, beq, nx * N * sizeof(double));
    memcpy(u_new, beq, nx * N * sizeof(double));

    // for (size_t i = 0; i < m; i++) printf("l_new[%d]=%.3f\n", i + 1, l_new[i]);
    // for (size_t i = 0; i < m; i++) printf("u_new[%d]=%.3f\n", i + 1, u_new[i]);
    osqp_update_bounds(work, l_new, u_new);

}

int main(int argc, char** argv) {
    // 1. Simulation parameters
    int T_sim = 100;
    if (argc > 1) {
        T_sim = atoi(argv[1]);
    }
    
    int max_iter0 = 100;
    int max_iter = max_iter0;
    double current_x[nx];
    memcpy(current_x, x0, sizeof(x0)); // Start from initial state defined in MPC_problem.h
    osqp_update_max_iter(&workspace, max_iter);

    // 2. OSQP already maintains workspace for warm start by default

    // 3. Print CSV Header
    #if PRINT_DEBUG
    printf("step");
    for (size_t i = 0; i < nx; i++) printf(",x%zu", i + 1);
    for (size_t i = 0; i < nu; i++) printf(",u%zu", i + 1);
    printf(",runtime_ms,J\n");
    #endif

    double total_runtime = 0.0;
    double cumulative_cost = 0.0;

    // 4. Simulation Loop
    for (int t = 0; t < T_sim; t++) {
        double u0[nu];

        // --- Step A: Get Control from OSQP ---
         // warm start in the first step, without clocking
        if (t == 0) {
            max_iter = 100000; // Allow more iterations in the first step for better convergence
        } else{
            max_iter = max_iter0;
        }

        max_iter = max_iter0;
        osqp_update_max_iter(&workspace, max_iter);

        /* 
         * Note: To properly update the initial state in OSQP, one should update 
         * the lower and upper bounds (l and u) at the indices corresponding 
         * to the initial state constraint. 
         */
        update_x0(&workspace, current_x);

        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
       
        osqp_solve(&workspace);
        
        clock_gettime(CLOCK_MONOTONIC, &end);
        double runtime_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
        if (t > 0){
            total_runtime += runtime_ms;
        }

        // Extract the first control input (assuming it starts at index 0 of the solution variables)
        for (size_t i = 0; i < nu; i++) {
            u0[i] = workspace.solution->x[i + N*nx];
        }

        double total_cost = 0.0;
        evaluate_cost_J(current_x, u0, &total_cost);
        cumulative_cost += total_cost;

        // --- Step B: Logging ---
        #if PRINT_DEBUG
        printf("%d", t);
        for (size_t i = 0; i < nx; i++) printf(",%.4f", current_x[i]);
        for (size_t i = 0; i < nu; i++) printf(",%.4f", u0[i]);
        printf(",%.4f,%.4f\n", runtime_ms, total_cost);
        #endif

        // --- Step C: Plant Simulation (Discrete Dynamics) ---
        // x(t+1) = A * x(t) + B * u(t)
        double x_plus[nx];
        pempc_aAx('n', 1.0, A, current_x, nx, nx, x_plus);
        pempc_aAxpy('n', 1.0, B, u0, nx, nu, x_plus);
        memcpy(current_x, x_plus, nx * sizeof(double));
    }

    // 5. Print Runtime Summary
    printf("Total runtime over %d steps: %.4f ms\n", T_sim, total_runtime);
    printf("Average runtime per step: %.4f ms\n", total_runtime / T_sim);
    printf("Cumulative cost: %.4f\n", cumulative_cost);

    return 0;
}
