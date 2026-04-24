
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <blasfeo_d_aux_ext_dep.h>

#include <hpipm_d_ocp_qp_ipm.h>
#include <hpipm_d_ocp_qp_dim.h>
#include <hpipm_d_ocp_qp.h>
#include <hpipm_d_ocp_qp_sol.h>
#include <hpipm_d_ocp_qp_utils.h>
#include <hpipm_timing.h>

#include "MPC_problem.h"


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

int n_box_u_cons = 0; // all input constrained
int n_box_x_cons = 0; // all states constrained
int box_u_idxs[nu];
int box_x_idxs[nx];
// the box constraint bounds without the Inf entries
double box_u_lb[nu];
double box_u_ub[nu];
double box_x_lb[nx];
double box_x_ub[nx];

void setup_box_constraint_dim_idx(){
    // set dimensions
    n_box_u_cons = 0; // all input constrained
    n_box_x_cons = 0; // all states constrained
    // box_u_idxs[nu];
    // box_x_idxs[nx];
    for (int ii = 0; ii < nu; ii++){
        if (umin[ii] != Inf && umax[ii] != Inf && umin[ii] != -Inf && umax[ii] != -Inf){
            box_u_idxs[n_box_u_cons] = ii;
            box_u_lb[n_box_u_cons] = umin[ii];
            box_u_ub[n_box_u_cons] = umax[ii];
            n_box_u_cons++;
        }
    }
    for (int ii = 0; ii < nx; ii++){
        if (xmin[ii] != Inf && xmax[ii] != Inf && xmin[ii] != -Inf && xmax[ii] != -Inf){
            box_x_idxs[n_box_x_cons] = ii;
            box_x_lb[n_box_x_cons] = xmin[ii];
            box_x_ub[n_box_x_cons] = xmax[ii];
            n_box_x_cons++;
        }
    }

}

// HPIPM is pretty complicated.
// We basically follow the routine from example_d_ocp_qp.c

/** 
* Set dimension data for the OCP QP; We take the getting_started_data.c & mass_spring_qp_data.c data as example
*/
void set_dimention_data(struct d_ocp_qp_dim* dim){
    setup_box_constraint_dim_idx();

    // number of input
    static int nnu[N+1];
    // number of states
    static int nnx[N+1];
    // number of input box constraints
    static int nnbu[N+1];
    // number of states box constraints
    static int nnbx[N+1];
    // number of general constraints
    static int nng[N+1];
    // number of slacks of soft constraints
    static int nns[N+1];
    // number of input box constraints considered as equality
    static int nnbue[N+1];
    // number of states box constraints considered as equality
    static int nnbxe[N+1];
    // number of general constraints considered as equality
    static int nnge[N+1];

    for (int ii=0; ii<=N; ii++){
        nnx[ii] = nx;
        nnu[ii] = nu;

        // number of state box constraints
        if (ii == 0) {
            nnbx[ii] = nx;
        }else{
            nnbx[ii] = n_box_x_cons;
        }

        // number of input box constraints
        nnbu[ii] = (ii < N) ? n_box_u_cons : 0;

        nng[ii]  = 0;
        nns[ii]  = 0;
        nnbue[ii]= 0;

        // x0 constraint is equality
        nnbxe[ii]= (ii == 0) ? (nx) : 0;
        nnge[ii] = 0;
    }

	// unified setter
	d_ocp_qp_dim_set_all(nnx, nnu, nnbx, nnbu, nng, nns, dim);

}


/**
* Set QP data for the OCP QP; We take the getting_started_data.c & mass_spring_qp_data.c data as example
We can also look at test_d_ocp.c for an example
*/
void set_qp_data(struct d_ocp_qp* qp){


    int idxbx_all[nx];
    for (int ii=0; ii<nx; ii++){
        idxbx_all[ii] = ii;
    }
    int idxbu_all[nu];
    for (int ii=0; ii<nu; ii++){
        idxbu_all[ii] = ii;
    }

    for (int ii = 0; ii < N; ii++) {
        // printf("Setting QP data i=%d\n", ii);
        d_ocp_qp_set_A(ii, A, qp);
        d_ocp_qp_set_B(ii, B, qp);
        // d_ocp_qp_set_b(ii, b_data, qp);
        d_ocp_qp_set_Q(ii, Q, qp);
        // d_ocp_qp_set_q(ii, q_data, qp);
        d_ocp_qp_set_R(ii, R, qp);

        if (ii == 0){
            d_ocp_qp_set_idxbx(ii, idxbx_all, qp);
        } else{
            d_ocp_qp_set_idxbx(ii, box_x_idxs, qp);
        }
        d_ocp_qp_set_idxbu(ii, box_u_idxs, qp);

        if (ii > 0) {
            d_ocp_qp_set_lbx(ii, box_x_lb, qp);
            d_ocp_qp_set_ubx(ii, box_x_ub, qp);
        }

        d_ocp_qp_set_lbu(ii, box_u_lb, qp);
        d_ocp_qp_set_ubu(ii, box_u_ub, qp);
    }

	// d_ocp_qp_set_all(
 //        hA, hB, hb, hQ, hS, hR, hq, hr, // The System dynamics & Cost
 //        hidxbx, hlbx, hubx, hidxbu, hlbu, hubu,
 //        hC, hD, hlg, hug, hZl, hZu, hzl, hzu, hidxs, hidxs_rev, hlls, hlus, qp);
}

void update_x0(struct d_ocp_qp* qp, double* x0_){
    // static int idxe0[nx];
    // for (int ii=0; ii<nx; ii++){
    //     idxe0[ii] = ii;
    // }
    // d_ocp_qp_set_idxe(0, idxe0, qp);

    d_ocp_qp_set_lbx(0, x0_, qp);
    d_ocp_qp_set_ubx(0, x0_, qp);
}

/***
* We Follow the example_d_ocp_qp.c file as a template
*/
int main(int argc, char **argv){

    int T_sim = 100;
    if (argc > 1) {
        T_sim = atoi(argv[1]);
    }

    /*********ocp qp dim
    ************************************************/

	hpipm_size_t dim_size = d_ocp_qp_dim_memsize(N);
	void *dim_mem = malloc(dim_size);

	struct d_ocp_qp_dim dim;
	d_ocp_qp_dim_create(N, &dim, dim_mem);

    set_dimention_data(&dim);

	// d_ocp_qp_dim_print(&dim);
    printf("OCP QP dimension set.\n");

    hpipm_size_t qp_size = d_ocp_qp_memsize(&dim);
	void *qp_mem = malloc(qp_size);

	struct d_ocp_qp qp;
	d_ocp_qp_create(&dim, &qp, qp_mem);
    set_qp_data(&qp);

    printf("OCP QP data set.\n");

    update_x0(&qp, x0);
	// d_ocp_qp_print(&dim, &qp);

    /******** ocp qp sol ********************************/

	hpipm_size_t qp_sol_size = d_ocp_qp_sol_memsize(&dim);
	void *qp_sol_mem = malloc(qp_sol_size);

	struct d_ocp_qp_sol qp_sol;
	d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    /********* ipm arg
    ************************************************/

	hpipm_size_t ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
	void *ipm_arg_mem = malloc(ipm_arg_size);

	struct d_ocp_qp_ipm_arg arg;
    int mode = 0; // copied from getting_started_data.c
	d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);
	d_ocp_qp_ipm_arg_set_default(mode, &arg);

    int warm_start = 1;
	d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);

    /******** ipm workspace
    ************************************************/

	hpipm_size_t ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
	void *ipm_mem = malloc(ipm_size);

	struct d_ocp_qp_ipm_ws workspace;
	d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    /****** ipm solver
    ************************************************/


    int hpipm_status;

    double current_x[nx];
    memcpy(current_x, x0, sizeof(x0)); // Start from initial state defined in MPC_problem.h

    // 3. Print CSV Header
    #if PRINT_DEBUG
    printf("step");
    for (size_t i = 0; i < nx; i++) printf(",x%zu", i + 1);
    for (size_t i = 0; i < nu; i++) printf(",u%zu", i + 1);
    printf(",runtime_ms,J\n");
    #endif


    // 4. Simulation Loop
    double total_runtime = 0.0;
    double cumulative_cost = 0.0;
    for (int t = 0; t < T_sim; t++){
        struct timespec start, end;
        update_x0(&qp, current_x);

        clock_gettime(CLOCK_MONOTONIC, &start);
	    d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
        clock_gettime(CLOCK_MONOTONIC, &end);

        double runtime_ms = (end.tv_sec - start.tv_sec) * 1000.0 + (end.tv_nsec - start.tv_nsec) / 1000000.0;
        if (t > 0){
            total_runtime += runtime_ms;
        }

	    d_ocp_qp_ipm_get_status(&workspace, &hpipm_status);
        // printf("\nHPIPM returned with flag %i.\n", hpipm_status);
        if(hpipm_status == 0){
            // printf("\n -> QP solved!\n");
        }
	    else if(hpipm_status==1){printf("\n -> Solver failed! Maximum number of iterations reached\n");}
	    else if(hpipm_status==2){printf("\n -> Solver failed! Minimum step lenght reached\n");}
	    else if(hpipm_status==3){printf("\n -> Solver failed! NaN in computations\n");}
	    else{printf("\n -> Solver failed! Unknown return flag\n");}
        // printf("\nAverage solution time over %i runs: %e [s]\n", nrep, time_ipm);
	    // printf("\n\n");
	    // printf("\nu = \n");

        double u0[nu];
	    d_ocp_qp_sol_get_u(0, &qp_sol, u0);

        double total_cost = 0.0;
        evaluate_cost_J(current_x, u0, &total_cost);
        cumulative_cost += total_cost;

        // --- Step B: Logging ---
        #if PRINT_DEBUG
        printf("%d", t);
        for (size_t i = 0; i < nx; i++) printf(",%.4f", current_x[i]);
        for (size_t i = 0; i < nu; i++) printf(",%.4f", u0[i]);
        printf(",%.4f,%4f\n", runtime_ms,total_cost);
        #endif

        // --- Step C: Plant Simulation (Forward Euler / Discrete Dynamics) ---
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


    /******** free memory and return ************************************************/
    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

	
    return 0;
}