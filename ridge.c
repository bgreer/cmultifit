#include "header.h"
#include "math.h"

/* Given a mode wavenumber and frequency, determing which ridge it is (f-mode, p1, p2, etc) */
int detect_ridge (double k, double nu)
{
	int ii, closest;
	double distance, temp;
	/* coefficients of polynomials in log-log space for each ridge */
	int known = 11;
	double c[] = {7.87, 0.48, 0.00, /* f-mode */
				  8.08, 0.47, 0.06, /* p1 */
				  8.26, 0.46, 0.05, /* p2 */
				  8.41, 0.44, 0.04, 
				  8.55, 0.45, 0.04, 
				  8.67, 0.47, 0.04, 
				  8.77, 0.47, 0.04, 
				  8.85, 0.47, 0.04, 
				  8.92, 0.47, 0.04, 
				  8.99, 0.47, 0.04, 
				  9.05, 0.47, 0.04  /* p10 */
				  };

	distance = 1e9;
	for (ii=0; ii<known; ii++)
	{
		if ((temp=distance_to_ridge(k, nu, c[ii*3], c[ii*3+1], c[ii*3+2])) < distance)
		{
			distance = temp;
			closest = ii;
		}
	}
	return closest;
}

double distance_to_ridge (double k, double nu, double a, double b, double c)
{
	double poly, logk;

	logk = log(k);
	poly = a + b*logk + c*logk*logk;

	return fabs(poly-log(nu));
}
