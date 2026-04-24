#include <stdio.h>
#include <string.h>
#include <time.h>
#include "RRLB_MPC.h"
#include "MPC_problem.h"
#include <stdlib.h> 

#define PRINT_DEBUG 1

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

/**
 * @brief Phase 3: Closed-loop simulation runner.
 * This program simulates the plant dynamics using the C-implemented MPC controller.
 * It outputs the state and control trajectories in CSV format for comparison with MATLAB.
 */
int main(int argc, char** argv) {
    // 1. Simulation parameters
    int T_sim = 100;
    if (argc > 1) {
        T_sim = atoi(argv[1]);
    }
    const size_t max_iter0 = 5;
    size_t max_iter = max_iter0;
    const double tol = 1e-4;
    
    double current_x[nx];
    memcpy(current_x, x0, sizeof(x0)); // Start from initial state defined in MPC_problem.h

    // 2. Initialize stacked variables for warm start
    double z_stacked[N * (nx + nu)];
    double lam_stacked[N * nx];
    memset(z_stacked, 0, sizeof(z_stacked));
    memset(lam_stacked, 0, sizeof(lam_stacked));

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
        double z_next[N * (nx + nu)];
        double lam_next[N * nx];
        // warm start
        if (t == 0) {
            max_iter = 50000;
        }else{
            max_iter = max_iter0;
        }

        // --- Step A: Get Control from C-MPC ---
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);
        mpc_get_control(current_x, max_iter, tol, z_stacked, lam_stacked, z_next, lam_next, u0);
        clock_gettime(CLOCK_MONOTONIC, &end);
        double runtime_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;

        double total_cost = 0.0;
        evaluate_cost_J(current_x, u0, &total_cost);
        cumulative_cost += total_cost;

        if (t > 0){
            total_runtime += runtime_ms;
        }

        // --- Step B: Logging ---
        #if PRINT_DEBUG
        printf("%d", t);
        for (size_t i = 0; i < nx; i++) printf(",%.4f", current_x[i]);
        for (size_t i = 0; i < nu; i++) printf(",%.4f", u0[i]);
        printf(",%.4f,%.4f\n", runtime_ms, total_cost);
        #endif

        // --- Step C: Plant Simulation (Forward Euler / Discrete Dynamics) ---
        // x(t+1) = A * x(t) + B * u(t)
        double x_plus[nx];
        pempc_aAx('n', 1.0, A, current_x, nx, nx, x_plus);
        pempc_aAxpy('n', 1.0, B, u0, nx, nu, x_plus);
        memcpy(current_x, x_plus, nx * sizeof(double));

        // --- Step D: Warm Start Update ---
        memcpy(z_stacked, z_next, sizeof(z_stacked));
        memcpy(lam_stacked, lam_next, sizeof(lam_stacked));
    }

    // 5. Print Runtime Summary
    printf("Total runtime over %d steps: %.4f ms\n", T_sim, total_runtime);
    printf("Average runtime per step: %.4f ms\n", total_runtime / T_sim);
    printf("Cumulative cost: %.4f\n", cumulative_cost);

    return 0;
}
