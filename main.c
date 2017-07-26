//////////////////////////////////////////////////////////////////////////////////////
//                                                                             		//
//  	Test case for linear system solvers			                               	//           
//                                                                             		//  
//  	Developer:  Jiho Yang (MEng)                                               	//
//              	M.Sc. candidate, Computational Science & Engineering           	//
//              	Technische Universitat Munchen                                 	//
//                                                                             		// 
//  	Final update date: 21/07/2017                                              	//
//                                                                             		//
//////////////////////////////////////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "helper.h"
#include "linear_solvers.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))

int main(int argn, char **args){   

    int itermax = 500;      
    int i, j;
	int n_rows = 3;
	int n_cols = 3;
	int dim = MAX(n_rows, n_cols);

    double **A;
    double *u;
	double *b;         

	double eps = 0.0001;
	double weight = 0.6;

	// Allocate memory
    A = matrix(0, n_rows-1, 0, n_cols-1);
	b = malloc((dim) * sizeof(double));
    u = malloc((dim) * sizeof(double));

	// Initialise matrix A
    init_matrix(A, 0, n_rows-1, 0, n_cols-1, 0);

	// Initialise solution and vector b
    for (i = 0; i < n_rows; i++){
        u[i] = 1;
		b[i] = 1;
		A[i][i] = 3;
    }

	// Linear Solver
	//weighted_jacobi(A, u, b, itermax, eps, n_rows, n_cols, dim, weight);
	sor(A, u, b, itermax, eps, n_rows, n_cols, dim, weight);
	//multigrid(A, u, b, itermax, eps, n_rows, n_cols, dim, cycle);
	//conjugate_gradient(A, u, b, itermax, eps, n_rows, n_cols, dim);

	// Print results
	printf("-----------Computed solution-----------\n\n");
	for (i = 0; i < dim; i++){
		printf("u[%i] = %f\n", i,u[i]);
	}
	printf("\n");

	// Print correct solution
	printf("-----------Backslash (matlab) solution-----------\n\n");
	printf("0.3333\n");
	printf("0.3333\n");
	printf("0.3333\n");
	printf("\n");

    // Free memory allocation
    free_matrix(A, 0, n_rows-1, 0, n_cols-1);
    free(u);
	free(b);

	printf("End of the program\n");

    return 0;
}

