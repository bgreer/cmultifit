#include "header.h"
#include <math.h>

int funk (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	int istw, iendw, num, iw, itht, ik;
	double akt, w1, tht, den, x2, err, back, back2, cost, sint;
	
	sub = (struct kslice*) private;
	istw = sub->start;
	iendw = sub->end;
	akt = sub->k;

	/* Add background first */
	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		back = p[n-3]*p[n-1]/((w1-p[n-2])*(w1-p[n-2]) + 0.25*p[n-1]*p[n-1]);
		back2 = 1.0/(1.0+pow(p[n-7]*w1, p[n-6]));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = TWOPI*itht/sub->ntheta;
			deviates[num] = back + back2*p[n-8]*(1.0 + p[n-5]*cos(2.0*(tht-p[n-4])));
			num++;
		}
	}

	/* Add each peak individually */
	for (ik=0; ik<(n-NBACK)/NPEAK; ik++)
	{
		num = 0;
		/*  pre-compute theta variation */
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = TWOPI*itht/sub->ntheta;
			thtarr[itht] = akt*(p[ik*NPEAK+3]*cos(tht)+p[ik*NPEAK+4]*sin(tht))/TWOPI - p[ik*NPEAK];
			thtpow[itht] = 0.5*p[ik*NPEAK+2]*p[ik*NPEAK+1]*(1.0+p[ik*NPEAK+5]*cos(2.0*(tht-p[ik*NPEAK+6])));
		}
		/* Loop over each nu, then theta */
		for (iw=istw; iw<=iendw; iw++)
		{
			w1 = iw*sub->delta_nu;
			for (itht=0; itht<sub->ntheta; itht++)
			{
				den = w1 + thtarr[itht];
				den = den*den + p[ik*NPEAK+2]*p[ik*NPEAK+2]/4.0;
				x2 = thtpow[itht]/den;
				deviates[num] += x2;
				num++;
			}
		}
	}

	/* Second step: subtract from data, weight appropriately */
	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		for (itht=0; itht<sub->ntheta; itht++)
		{
			if (sub->data[iw-istw][itht] > 0.0 && !isnan(sub->data[iw-istw][itht]))
			{
				/* Weight Chi appropriately */
				switch (sub->par->chiweight)
				{
					case WEIGHT_NOISE:
						deviates[num] = (deviates[num] - sub->data[iw-istw][itht]);
						deviates[num] /= sub->noise[iw-istw][itht];
						break;
					case WEIGHT_MAXL:
						deviates[num] = sub->data[iw-istw][itht]/deviates[num];
						deviates[num] -= log(deviates[num]);
						deviates[num] = sqrt(deviates[num]);
						break;
				}
			} else {
				deviates[num] = 0.0;
			}
			num++;
		}
	}
	return EXIT_SUCCESS;
}
