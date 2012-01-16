#include "header.h"
#include <math.h>

double ml_funk (const gsl_vector *v, void *params)
{
	/* v contains starting parameters */
	/* params is a pointer to a struct kslice */
	struct kslice *sub;
	double akt, w1, tht, den, x2, back, back2, backw, freq;
	int ii, iw, istw, iendw, n, itht, ij;
	double model, lh, lh2;
	sub = (struct kslice*) params;

	lh = 0.0;
	lh2 = 0.0;
	istw = sub->start;
	iendw = sub->end;
	akt  = sub->k;
	n = sub->n;

	for (iw = istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		back = 0.5*gvg(v,n*NPEAK+5)*gvg(v,n*NPEAK+7) / 
			((w1-gvg(v,n*NPEAK+6))*(w1-gvg(v,n*NPEAK+6)) + 0.25*gvg(v,n*NPEAK+7)*gvg(v,n*NPEAK+7));
		back2 = 1.0/(1.0+pow(w1/gvg(v,n*NPEAK+1), gvg(v,n*NPEAK+2)));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			if (sub->data[iw-istw][itht] <= 0.0)
			{
				printf("ERROR\n");
				exit(0);
			}
			model = back;
			tht = TWOPI*itht/sub->ntheta;
			model += back2*gvg(v,n*NPEAK)*(1.0+gvg(v,n*NPEAK+3)*cos(2.0*(tht-gvg(v,n*NPEAK+4))));
			/* create and sum each peak */
			for (ii=0; ii<n; ii++)
			{
				den = w1+akt*(gvg(v,ii*NPEAK+3)*cos(tht)+gvg(v,ii*NPEAK+4)*sin(tht))/TWOPI
					- gvg(v,ii*NPEAK);
				den = den*den + 0.25*gvg(v,ii*NPEAK+2)*gvg(v,ii*NPEAK+2);
				model += (0.5*gvg(v,ii*NPEAK+2)*gvg(v,ii*NPEAK+1)*
					(1.0+gvg(v,ii*NPEAK+5)*cos(2.0*(tht-gvg(v,ii*NPEAK+6)))))/den;
			}

			lh += log(model/sub->data[iw-istw][itht]);
			lh2 += sub->data[iw-istw][itht]/model;
		}
	}

	/* Various shady penalties */
/*	for (ii=0; ii<n; ii++)
		if (gvg(v,ii*NPEAK) > sub->end*sub->delta_nu) lh2 += 20.*(gvg(v,ii*NPEAK)-sub->end*sub->delta_nu);
	for (ii=1; ii<n; ii++)
		if (gvg(v,ii*NPEAK) < gvg(v,(ii-1)*NPEAK)) lh2 += 100.*(gvg(v,(ii-1)*NPEAK)-gvg(v,ii*NPEAK));
	if (gvg(v,n*NPEAK+7) < 500.) lh2 += 100.*(gvg(v,n*NPEAK+7)-500.)*(gvg(v,n*NPEAK+7)-500.);
	for (ii=0; ii<n; ii++)
		if (gvg(v,ii*NPEAK+1) < 0.0) lh2 += 1000.*gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+1);
*/
	return lh + lh2;
}

void ml_dfunk (const gsl_vector *v, void *params, gsl_vector *g)
{
	struct kslice *sub;
	double akt, w1, tht, den, den1, x2, back, back2, backw, freq;
	int ii, iw, istw, iendw, n, itht, ij;
	double model, lor, amp, co, corr;
	double *deriv;
	sub = (struct kslice*) params;

	istw = sub->start;
	iendw = sub->end;
	akt  = sub->k;
	n = sub->n;
	deriv = calloc(n*NPEAK+NBACK, sizeof(double));

	/* Go through each peak first */
	for (iw = istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		back = gvg(v,n*NPEAK)/(1.+pow(w1/gvg(v,n*NPEAK+1), gvg(v,n*NPEAK+2)));
		back2 = 0.5*gvg(v,n*NPEAK+5)*gvg(v,n*NPEAK+7) 
			/ ((w1-gvg(v,n*NPEAK+6))*(w1-gvg(v,n*NPEAK+6)) + 0.25*gvg(v,n*NPEAK+7)*gvg(v,n*NPEAK+7));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = TWOPI*itht/sub->ntheta;
			backw = back*(1.+gvg(v,n*NPEAK+3)*cos(2.*(tht-gvg(v,n*NPEAK+4))));
			/* Create model */
			model = back2+backw;
			for (ii=0; ii<n; ii++)
			{
				den = w1 - gvg(v,ii*NPEAK)
					+ akt*(gvg(v,ii*NPEAK+3)*cos(tht)+gvg(v,ii*NPEAK+4)*sin(tht))/TWOPI;
				den = den*den + 0.25*gvg(v,ii*NPEAK+2)*gvg(v,ii*NPEAK+2);
				model += 0.5*gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+2) 
					*(1+gvg(v,ii*NPEAK+5)*cos(2.*(tht-gvg(v,ii*NPEAK+6))))/den;
			}
			corr = (1.-sub->data[iw-istw][itht]/model)/model;

			/* Calculate each derivative */
			for (ii=0; ii<n; ii++)
			{
				co = cos(2.*(tht-gvg(v,ii*NPEAK+6)));
				den = w1 - gvg(v,ii*NPEAK)
					+ akt*(gvg(v,ii*NPEAK+3)*cos(tht)+gvg(v,ii*NPEAK+4)*sin(tht))/TWOPI;
				lor = 0.5*gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+2)*(1.+gvg(v,ii*NPEAK+5)*co);
				lor /= den*den + 0.25*gvg(v,ii*NPEAK+2)*gvg(v,ii*NPEAK+2);

				/* Central frequency */
				deriv[ii*NPEAK+0] += corr*lor*4.*lor*den 
					/(gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+2)*(1.+gvg(v,ii*NPEAK+5)*co));
				/* Amplitude */
				deriv[ii*NPEAK+1] += corr*lor/gvg(v,ii*NPEAK+1);
				/* Width */
				deriv[ii*NPEAK+2] += corr*lor*(1.-lor*gvg(v,ii*NPEAK+2)/(1.+gvg(v,ii*NPEAK+5)*co) 
					/gvg(v,ii*NPEAK+1))/gvg(v,ii*NPEAK+2);
				/* Ux */
				deriv[ii*NPEAK+3] += -corr*lor*lor*den*(akt*cos(tht)/PI)*2. 
					/(gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+2)*(1.+gvg(v,ii*NPEAK+5)*co));
				/* Uy */
				deriv[ii*NPEAK+4] += -corr*lor*lor*den*(akt*sin(tht)/PI)*2.
					/(gvg(v,ii*NPEAK+1)*gvg(v,ii*NPEAK+2)*(1.+gvg(v,ii*NPEAK+5)*co));
				/* Anisotropy fraction */
				deriv[ii*NPEAK+5] += corr*lor*co/(1.+gvg(v,ii*NPEAK+5)*co);
				/* Anisotropy direction */
				deriv[ii*NPEAK+6] += corr*lor*2.*gvg(v,ii*NPEAK+5)*sin(2.*(tht-gvg(v,ii*NPEAK+6))) 
					/(1.+gvg(v,ii*NPEAK+5)*co);
			}
			
			/* Then the background terms */

			/* Power law amplitude */
			deriv[n*NPEAK+0] += corr*backw/gvg(v,n*NPEAK+0);
			/* Power law cut-off frequency */
			deriv[n*NPEAK+1] += corr*backw*back*gvg(v,n*NPEAK+2)*pow(w1/gvg(v,n*NPEAK+1),gvg(v,n*NPEAK+2)) 
				/ (gvg(v,n*NPEAK+0)*gvg(v,n*NPEAK+1));
			/* Power law index */
			if (w1 > 0.0)
				deriv[n*NPEAK+2] += -corr*backw*back*log(w1/gvg(v,n*NPEAK+1)) 
					* pow(w1/gvg(v,n*NPEAK+1), gvg(v,n*NPEAK+2)) / gvg(v,n*NPEAK);
			/* Power law anisotropy fraction */
			deriv[n*NPEAK+3] += corr*back*cos(2.*(tht-gvg(v,n*NPEAK+4)));
			/* Power law anisotropy direction */
			deriv[n*NPEAK+4] += corr*back*2.0*gvg(v,n*NPEAK+3)*sin(2.*(tht-gvg(v,n*NPEAK+4)));

			/* Lorentzian amplitude */
			deriv[n*NPEAK+5] += corr*back2/gvg(v,n*NPEAK+5);
			/* Lorentizian central frequency */
			deriv[n*NPEAK+6] += corr*back2*back2*4.*(w1-gvg(v,n*NPEAK+6))/(gvg(v,n*NPEAK+5)*gvg(v,n*NPEAK+7));
			/* Lorentzian width */
			deriv[n*NPEAK+7] += corr*back2*(1.-back2*gvg(v,n*NPEAK+7)/gvg(v,n*NPEAK+5))/gvg(v,n*NPEAK+7);
		}
	}


	for (ii=0; ii<n*NPEAK+NBACK; ii++)
		gsl_vector_set(g,ii,deriv[ii]);
	free(deriv);
}

void ml_fdfunk (const gsl_vector *v, void *params, double *f, gsl_vector *g)
{
	*f = ml_funk(v, params);
	ml_dfunk(v, params, g);
}

int funk (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	double akt, w1, tht, den, x2, back, back2, backw, freq, co, corr, lor;
	int ii, iw, istw, iendw, itht, ij, nr, num;
	double model, lh, lh2;
	sub = (struct kslice*) private;

	lh = 0.0;
	lh2 = 0.0;
	istw = sub->start;
	iendw = sub->end;
	akt  = sub->k;
	nr = sub->n;

	num = 0;
	for (iw = istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		back = p[nr*NPEAK]/(1.+pow(w1/p[nr*NPEAK+1], p[nr*NPEAK+2]));
		back2 = 0.5*p[nr*NPEAK+5]*p[nr*NPEAK+7] 
			/ ((w1-p[nr*NPEAK+6])*(w1-p[nr*NPEAK+6]) + 0.25*p[nr*NPEAK+7]*p[nr*NPEAK+7]);
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = TWOPI*itht/sub->ntheta;
			backw = back*(1.+p[nr*NPEAK+3]*cos(2.*(tht-p[nr*NPEAK+4])));
			/* Create model */
			model = back2+backw;
			for (ii=0; ii<nr; ii++)
			{
				den = w1 - p[ii*NPEAK]
					+ akt*(p[ii*NPEAK+3]*cos(tht)+p[ii*NPEAK+4]*sin(tht))/TWOPI;
				den = den*den + 0.25*p[ii*NPEAK+2]*p[ii*NPEAK+2];
				model += 0.5*p[ii*NPEAK+1]*p[ii*NPEAK+2] 
					*(1+p[ii*NPEAK+5]*cos(2.*(tht-p[ii*NPEAK+6])))/den;
			}
			lh = model / sub->data[iw-istw][itht];
			if (sub->data[iw-istw][itht] > 0.0)
			{
				deviates[num] = 0.5*(model-sub->data[iw-istw][itht])/sqrt(model*sub->data[iw-istw][itht]);
				corr = (1.-sub->data[iw-istw][itht]/model)/model;
			}
			else
			{
				deviates[num] = 0.0;
				corr = 0.0;
			}
			if (derivs)
			{
				/* Calculate each derivative */
				for (ii=0; ii<nr; ii++)
				{
					co = cos(2.*(tht-p[ii*NPEAK+6]));
					den = w1 - p[ii*NPEAK]
						+ akt*(p[ii*NPEAK+3]*cos(tht)+p[ii*NPEAK+4]*sin(tht))/TWOPI;
					lor = 0.5*p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*co);
					lor /= den*den + 0.25*p[ii*NPEAK+2]*p[ii*NPEAK+2];
					/* Central frequency */
					if (derivs[ii*NPEAK+0])
					{
						derivs[ii*NPEAK+0][num] = corr*lor*4.*lor*den 
							/(p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*co));
					}
					/* Amplitude */
					if (derivs[ii*NPEAK+1])
						derivs[ii*NPEAK+1][num] = corr*lor/p[ii*NPEAK+1];
					/* Width */
					if (derivs[ii*NPEAK+2])
						derivs[ii*NPEAK+2][num] = corr*lor*(1.-lor*p[ii*NPEAK+2]/(1.+p[ii*NPEAK+5]*co) 
							/p[ii*NPEAK+1])/p[ii*NPEAK+2];
					/* Ux */
					if (derivs[ii*NPEAK+3])
						derivs[ii*NPEAK+3][num] = -corr*lor*lor*den*(akt*cos(tht)/PI)*2. 
							/(p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*co));
					/* Uy */
					if (derivs[ii*NPEAK+4])
						derivs[ii*NPEAK+4][num] = -corr*lor*lor*den*(akt*sin(tht)/PI)*2.
							/(p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*co));
					/* Anisotropy fraction */
					if (derivs[ii*NPEAK+5])
						derivs[ii*NPEAK+5][num] = corr*lor*co/(1.+p[ii*NPEAK+5]*co);
					/* Anisotropy direction */
					if (derivs[ii*NPEAK+6])
						derivs[ii*NPEAK+6][num] = corr*lor*2.*p[ii*NPEAK+5]*sin(2.*(tht-p[ii*NPEAK+6])) 
							/(1.+p[ii*NPEAK+5]*co);
				}
			
				/* Then the background terms */
				/* Power law amplitude */
				if (derivs[nr*NPEAK+0])
					derivs[nr*NPEAK+0][num] = corr*backw/p[nr*NPEAK+0];
				/* Power law cut-off frequency */
				if (derivs[nr*NPEAK+1])
					derivs[nr*NPEAK+1][num] = corr*backw*back*p[nr*NPEAK+2]*pow(w1/p[nr*NPEAK+1],p[nr*NPEAK+2]) 
						/ (p[nr*NPEAK+0]*p[nr*NPEAK+1]);
				/* Power law index */
				if (w1 > 0.0)
					if (derivs[nr*NPEAK+2])
						derivs[nr*NPEAK+2][num] = -corr*backw*back*log(w1/p[nr*NPEAK+1]) 
							* pow(w1/p[nr*NPEAK+1], p[nr*NPEAK+2]) / p[nr*NPEAK];
				/* Power law anisotropy fraction */
				if (derivs[nr*NPEAK+3])
					derivs[nr*NPEAK+3][num] = corr*back*cos(2.*(tht-p[nr*NPEAK+4]));
				/* Power law anisotropy direction */
				if (derivs[nr*NPEAK+4])
					derivs[nr*NPEAK+4][num] = corr*back*2.0*p[nr*NPEAK+3]*sin(2.*(tht-p[nr*NPEAK+4]));

				/* Lorentzian amplitude */
				if (derivs[nr*NPEAK+5])
					derivs[nr*NPEAK+5][num] = corr*back2/p[nr*NPEAK+5];
				/* Lorentizian central frequency */
				if (derivs[nr*NPEAK+6])
					derivs[nr*NPEAK+6][num] = corr*back2*back2*4.*(w1-p[nr*NPEAK+6])/(p[nr*NPEAK+5]*p[nr*NPEAK+7]);
				/* Lorentzian width */
				if (derivs[nr*NPEAK+7])
					derivs[nr*NPEAK+7][num] = corr*back2*(1.-back2*p[nr*NPEAK+7]/p[nr*NPEAK+5])/p[nr*NPEAK+7];

			}
			num++;
		}
	}


	return EXIT_SUCCESS;
}

