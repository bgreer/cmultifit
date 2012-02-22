#include "header.h"
#include <math.h>

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
					*(1.+p[ii*NPEAK+5]*cos(2.*(tht-p[ii*NPEAK+6])))/den;
			}
			lh = model / sub->data[iw-istw][itht];
			if (sub->data[iw-istw][itht] > 0.0)
			{
				/* likelihood = log(model/data) + data/model */
				/*deviates[num] = 0.5*(model-sub->data[iw-istw][itht])/sqrt(model*sub->data[iw-istw][itht]);*/
				deviates[num] = sqrt(log(model/sub->data[iw-istw][itht]) + sub->data[iw-istw][itht]/model);
				/*deviates[num] = 0.7071*(model-sub->data[iw-istw][itht])/model;*/
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
						derivs[ii*NPEAK+0][num] = corr*lor*4.*lor*den 
							/(p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*co));
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

