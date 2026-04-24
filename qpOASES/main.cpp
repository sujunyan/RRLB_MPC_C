#include <qpOASES.hpp>
#include <stdio.h>
#include <string.h>
#include <time.h>
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

void update_x0(double* beq, double* current_x) {
    double current_xk[nx];
    memcpy(current_xk, current_x, nx * sizeof(double));

    for (size_t k = 0; k < N; k++) {
        double next_xk[nx];
        pempc_aAx('n', 1.0, A, current_xk, nx, nx, next_xk);
        memcpy(current_xk, next_xk, nx * sizeof(double));

        for (size_t i = 0; i < nx; i++) {
            beq[k * nx + i] = current_xk[i];
        }
    }
}

void col_major_to_row_major(const double* input, double* output, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            output[i * cols + j] = input[j * rows + i];
        }
    }
}

int main(int argc, char** argv) {
    /*
     * @brief Phase 3: Closed-loop simulation runner for qpOASES.
     * This program simulates the plant dynamics using the qpOASES-implemented MPC controller.
     * It outputs the state and control trajectories in CSV format for comparison with MATLAB.
     */
    USING_NAMESPACE_QPOASES

    // 1. Simulation parameters
    int T_sim = 100;
    if (argc > 1) {
        T_sim = atoi(argv[1]);
    }

    double qp_H[QP_nz * QP_nz];
    double qp_A[QP_ncons * QP_nz];

    // Note: qpOASES uses row-major storage, while the data defined in MPC_problem.h is in column-major.
    col_major_to_row_major(QP_H, qp_H, QP_nz, QP_nz);
    col_major_to_row_major(QP_A, qp_A, QP_ncons, QP_nz);

    double * qp_g = QP_g;
    double * qp_l = QP_l;
    double * qp_u = QP_u;

    double current_x[nx];
    memcpy(current_x, x0, sizeof(x0)); // Start from initial state defined in MPC_problem.h

    QProblem prob(QP_nz, QP_ncons);

    Options options;
    options.printLevel = PL_LOW;
    prob.setOptions(options);

    // Initialize
    int_t nWSR = 100000;
    prob.init(qp_H, qp_g, qp_A, NULL, NULL, qp_l, qp_u, nWSR);

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

        // --- Step A: Get Control from qpOASES ---
        struct timespec start, end;
        clock_gettime(CLOCK_MONOTONIC, &start);

        // Update bounds for new x0
        double beq[nx * N];
        update_x0(beq, current_x);

        double lbA_new[QP_ncons];
        double ubA_new[QP_ncons];
        memcpy(lbA_new, qp_l, QP_ncons * sizeof(double));
        memcpy(ubA_new, qp_u, QP_ncons * sizeof(double));
        memcpy(lbA_new, beq, nx * N * sizeof(double));
        memcpy(ubA_new, beq, nx * N * sizeof(double));

        // Hotstart for updated bounds
        nWSR = 100;
        prob.hotstart(qp_g, NULL, NULL, lbA_new, ubA_new, nWSR);

        clock_gettime(CLOCK_MONOTONIC, &end);
        double runtime_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
        if (t > 0){
            total_runtime += runtime_ms;
        }

        // Extract the first control input (assuming it starts at index 0 of the solution variables)
        real_t zOpt[QP_nz];
        prob.getPrimalSolution(zOpt);
        for (size_t i = 0; i < nu; i++) {
            u0[i] = zOpt[i + N * nx];
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