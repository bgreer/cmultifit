#include "header.h"

/* TODO: write this */

/*
	Needs pointer to final parameter values
	Compute 2nd derivatives of 'chi-squared' wrt every variable
	 (plus cross derivatives)
	Invert this matrix, output as covariance.
*/
/*
double*** covariance (int numridges, double *p, double ***spec)
{
	double** A;
	int ii;

	A = malloc((NPEAK*numridges+NBACK)*sizeof(double*));
	for (ii=0; ii<(NPEAK*numridges+NBACK); ii++)
		A[ii] = malloc((NPEAK*numridges+NBACK)*sizeof(double));
	
	hessian (&A, p, spec);
	return &A;
}
*/
/* 
	A = empty square matrix
	p = final parameter values
	spec = spectrum
*/
void hessian (double ***A, double *p, double ***spec)
{
	
}
