//////////////////////////////////////////////////////////////////////////////
//																			//
//		Linear Solver Library for solving System of Linear Equations		//
//		Parallelised with OpenMP for shared memory systems					//
//																			//
//		Developer: Jiho Yang (MEng)											//
//				   M.Sc. candidate, Computational Science & Engineering 	//
//				   Technische Universitat Munchen							//
//																			//
//		Final Update: 22/07/2017											//
//																			//
//		TODO:																//
//			 1. Add Multigrid method										//
//			 2. Add Conjugate Gradient method								//
//			 3. Code optimisation using loop unrolling						//
//			 4. Optimise cache use											//
//			 5. Parallelise with OpenMP										//
//			 6. Add as template												//
//			 7. Consider classes on CPP										//
//			 8. Consider parallelisation with CUDA							//
//			 9. Add GMRES													//
//																			//
//////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linear_solvers.h"
#include "helper.h"
//#include <omp.h>

// Residual computation
double compute_residual(
					   	double **A,
					   	double *u,
					   	double *b,
						int n_rows,
						int n_cols,
						int dim
					   ){

		double Au = 0;
		double rloc = 0;
		double res = 0;
		int i, j;

		for (i = 0; i < n_rows; ++i){
			Au = 0;
			for (j = 0; j < n_cols; ++j){
				Au += A[i][j] * u[j];		
			}
			rloc += (b[i] - Au) * (b[i] - Au);
		}
		res = sqrt(rloc / dim);
		return res;
}

// Weighted Jacobi Relaxation
void weighted_jacobi(
					 double **A,
					 double *u,
					 double *b,
					 int itermax,
					 double eps,
					 int n_rows,
					 int n_cols,
					 int dim,
					 double weight
					){

	// Initialisation
	double res = 1;
	double Au = 0;
	int it = 0;
	int i,j;

	// Memory allocation for u_new
	double *u_new = malloc(dim * sizeof(double));

	// Main iteration
	while (it < itermax && res > eps){
		for (i = 0; i < n_rows; ++i){
			for (j = 0; j < i; ++j){
				Au += A[i][j] * u[j];
			}
			for (j = i+1; j < n_cols; ++j){
				Au += A[i][j] * u[j];
			}
			u_new[i] = weight * (b[i] - Au) / A[i][i] + (1 - weight) * u[i];
		}
		for (i = 0; i < n_rows; ++i){
			u[i] = u_new[i];
		}
		res = compute_residual(A, u, b, n_rows, n_cols, dim);
		it++;
	}
	printf("it = %i\n", it);
	free(u_new);
}

// Successive over-relaxation
void sor(
		 double **A,
		 double *u,
		 double *b,
		 int itermax,
		 double eps,
		 int n_rows,
		 int n_cols,
		 int dim,
		 double weight
		){

	// Initialisation
	double res = 1;
	double Au = 0;
	int it = 0;
	int i,j;

	// Main iteration
	while (it < itermax && res > eps){
		for (i = 0; i < n_rows; ++i){
			for (j = 0; j < i; ++j){
				Au += A[i][j] * u[j];
			}
			for (j = i+1; j < n_cols; ++j){
				Au += A[i][j] * u[j];
			}
			u[i] = weight * (b[i] - Au) / A[i][i] + (1 - weight) * u[i];
		}
		res = compute_residual(A, u, b, n_rows, n_cols, dim);
		it++;
	}
	printf("it = %i\n", it);
}

// Geometric Multigrid V cycle
void v_multigrid(
			   double **A,
			   double *u,
			   double *b,
			   int itermax,
			   double eps,
			   int n_rows,
			   int n_cols,
			   int dim,
			   int level_max
			  ){

}

// Geometric Multigrid W cycle
void w_multigrid(
			   double **A,
			   double *u,
			   double *b,
			   int itermax,
			   double eps,
			   int n_rows,
			   int n_cols,
			   int dim,
			   int level_max
			  ){

}

// Full Geometric Multigrid cycle (FMA)
void fma_multigrid(
			   double **A,
			   double *u,
			   double *b,
			   int itermax,
			   double eps,
			   int n_rows,
			   int n_cols,
			   int dim,
			   int level_max
			  ){

}

// Conjugate Gradient Method
void conjugate_gradient(
					    double **A,
					    double *u,
					    double *b,
					    int itermax,
						double eps,
					    int n_rows,
				   		int n_cols,
						int dim
					   ){

	double res = 1;
	double *r = malloc(dim * sizeof(double)); 
	int i,j;
	int it = 0;

	while (res > eps && it < itermax){
		double alpha = 0;
		// compute steepest descent
		for (i=0; i<n_rows; i++){
			double Au = 0;
			for (j=0; j<n_cols; j++){
				Au += A[i][j]*u[j];
			}
			r[i] = b[i] - Au;
		}
		// compute alpha
		for (i=0; i<n_rows; i++){
			double rA = 0;
			for (j=0; j<n_cols; j++){
				rA += r[j]*A[i][j];
			}
			alpha += (r[i]*r[i])/(rA*r[i]);
		}
		// update u
		for (i=0; i<n_rows; i++){
			u[i] = u[i] + alpha*r[i];
		}
		// compute residual
		res = compute_residual(A, u, b, n_rows, n_cols, dim);
		it++;
	}
	printf("it = %i\n", it);
}

// Generalized Minimal Residual Method (GMRES)
void gmres(
		   double **A,
		   double *u,
		   double *b,
		   int itermax,
		   double eps,
		   int n_rows,
		   int n_cols,
		   int dim
		  ){

}
