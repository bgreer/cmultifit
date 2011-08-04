#include "header.h"
#include "lapacke.h"
#include "lapacke_utils.h"

/* TODO: write this */

/*
	Needs a function pointer to the residual function
	and an array of best-fit parameters.
	Assume this array represents a minimum in parameter space.
	Numerically find the 2nd derivative of something

*/

void test ()
{
	lapack_int info, lda, n;
	char uplo;
	double *a;

	uplo = 'L';
	lda = 8;
	n = 4;

	a = (double*) malloc((n*lda)*sizeof(double));
	a[0] = 4.16000000000000010e+000;  /* a[0,0] */
    a[8] = 0.00000000000000000e+000;  /* a[0,1] */
    a[16] = 0.00000000000000000e+000;  /* a[0,2] */
    a[24] = 0.00000000000000000e+000;  /* a[0,3] */
    a[1] = -3.12000000000000010e+000;  /* a[1,0] */
    a[9] = 5.03000000000000020e+000;  /* a[1,1] */
    a[17] = 0.00000000000000000e+000;  /* a[1,2] */
    a[25] = 0.00000000000000000e+000;  /* a[1,3] */
    a[2] = 5.60000000000000050e-001;  /* a[2,0] */
    a[10] = -8.29999999999999960e-001;  /* a[2,1] */
    a[18] = 7.60000000000000010e-001;  /* a[2,2] */
    a[26] = 0.00000000000000000e+000;  /* a[2,3] */
    a[3] = -1.00000000000000010e-001;  /* a[3,0] */
    a[11] = 1.17999999999999990e+000;  /* a[3,1] */
    a[19] = 3.40000000000000020e-001;  /* a[3,2] */
    a[27] = 1.17999999999999990e+000;  /* a[3,3] */
	LAPACK_dpotrf(&uplo, &n, a, &lda, &info);


	printf("result? %d\n", info);
}
