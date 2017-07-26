#ifndef __LINEAR_SOLVERS_H_
#define __LINEAR_SOLVERS_H_	

double compute_residual(
						double **A,
						double *u,
						double *b,
						int n_rows,
						int n_cols,
						int dim
					   );

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
					);

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
		);

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
			    );

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
			    );

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
			    );

void conjugate_gradient(
						double **A,
						double *u,
						double *b,
						int itermax,
						double eps,
						int n_rows,
						int n_cols,
					    int dim
					   );

void gmres(
		   double **A,
		   double *u,
		   double *b,
		   int itermax,
		   double eps,
		   int n_rows,
		   int n_cols,
		   int dim
		  );

#endif
