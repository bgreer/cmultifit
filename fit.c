#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpfit.h"
#include "header.h"

int fit_peak (struct params* p, double *freq, double *amp, double *width, double*** pol, double delta_nu, double delta_k, int nnu, int ntheta, int k)
{
	int ii, ij, ik, mpreturn;
	double *param, *xerror;
	mp_par *bounds;
	mp_result *mpres2;
	mp_config *mpconf;
	struct kslice sub;

	param = (double*) malloc(4*sizeof(double));
	param[0] = *freq;
	param[1] = *amp;
	param[2] = *width;
	param[3] = (*amp)*0.1;

	mpres2 = (mp_result*) calloc(1,sizeof(mp_result));
	mpconf = (mp_config*) calloc(1, sizeof(mp_config));

	/* Load subsection data */
	sub.par = p;
	sub.start = (param[0] - 1.0*param[2])/delta_nu;
	sub.end = (param[0] + 1.0*param[2])/delta_nu;
	if (sub.start < 0) sub.start = 0;
	if (sub.end >= nnu) sub.end = nnu-1;
	sub.ntheta = ntheta;
	sub.delta_nu = delta_nu;
	sub.k = (k+1)*delta_k;
	sub.data = malloc((sub.end-sub.start+1)*sizeof(double*));
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		sub.data[ii-sub.start] = malloc(ntheta*sizeof(double));
		for (ik=0; ik<ntheta; ik++)
		{
			sub.data[ii-sub.start][ik] = pol[ii][k][ik];
		}
	}

	bounds = malloc(4*sizeof(mp_par));
	xerror = calloc(4, sizeof(double));

	for (ii=0; ii<4; ii++)
	{
		bounds[ii].fixed = 0;
		bounds[ii].side = 0;
		bounds[ii].deriv_debug = 0;
		bounds[ii].relstep = 0.0;
		bounds[ii].step = 0.0;
		bounds[ii].limited[0] = bounds[ii].limited[1] = 0;
	}
	bounds[0].limited[0] = bounds[0].limited[1] = 1;
	bounds[0].limits[0] = sub.start*delta_nu;
	bounds[0].limits[1] = sub.end*delta_nu;
	bounds[1].limited[0] = 1;
	bounds[1].limits[0] = 0.0;
	bounds[2].limited[0] = bounds[2].limited[1] = 1;
	bounds[2].limits[0] = param[2]/2.;
	bounds[2].limits[1] = param[2]*2.;
	bounds[3].limited[0] = 1;
	bounds[3].limited[1] = 0;
	bounds[3].limits[0] = 0.0;

	mpconf->ftol = 1e-8;
	mpconf->xtol = 1e-4;
	mpconf->gtol = 1e-8;
	mpconf->covtol = 1e-10;
	mpconf->maxiter = 200;
	mpconf->nofinitecheck = 1;

	mpres2->xerror = xerror;
	mpreturn = 1;
	mpreturn = mpfit(&funk_single, ntheta*(sub.end-sub.start+1), 4, 
					param, bounds, mpconf, &sub, mpres2);
	*freq = param[0];
	*amp = param[1]+param[3];/**20.;*param[0];*/
	*width = param[2];

	FILE *fp;
	double den;
	fp = fopen("debug2", "w");
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		for (ij=0; ij<ntheta; ij++)
		{
		den = ii*delta_nu - *freq;
		den = den*den + (*width)*(*width)/4.0;
		fprintf(fp, "%f\t%f\t%f\n", ii*delta_nu, param[3]+(*amp)*(*width)/(2.0*den), sub.data[ii-sub.start][ij]);
		}
	}
	fclose(fp);

	free(param);
	free(bounds);
	free(xerror);
	free(mpres2);
	free(mpconf);
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		free(sub.data[ii-sub.start]);
	}
	free(sub.data);

	return 0;
}

int funk_single (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	int istw, iendw, num, iw, itht, ik;
	double akt, w1, tht, den, x2, err;

	sub = (struct kslice*) private;
	istw = sub->start;
	iendw = sub->end;
	akt = sub->k;

	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = TWOPI*itht/sub->ntheta;
			deviates[num] = 0.0;
			if (sub->data[iw-istw][itht] > 0.0 && !isnan(sub->data[iw-istw][itht]))
			{
				den = w1 - p[0];
				den = den*den + p[2]*p[2]/4.0;
				x2 = p[1]*p[2]/(2.0*den)+p[3];
				deviates[num] = (x2 - sub->data[iw-istw][itht])/x2;
			}
			num++;
		}
	}
	return 0;
}

int fit_back (struct params* p, double* amp, double* cutoff, double* power, double*** pol, double delta_nu, int nnu, int ntheta, int k)
{
	int ii, ij, ik, mpreturn;
	double *param;
	mp_par *bounds;
	mp_result *mpres2;
	mp_config *mpconf;
	struct kslice sub;

	param = (double*) malloc(3*sizeof(double));
	
	/* Determine starting guesses */
	param[0] = 10.0;
	param[1] = 100.;
	param[2] = 2.5;

	mpres2 = (mp_result*) calloc(1,sizeof(mp_result));
	mpconf = (mp_config*) calloc(1, sizeof(mp_config));

	/* Load subsection data */
	sub.par = p;
	sub.start = (300.)/delta_nu;
	sub.end = (1500.)/delta_nu;
	if (sub.end >= nnu) sub.end = nnu-1;
	sub.ntheta = ntheta;
	sub.delta_nu = delta_nu;
	sub.k = 0.0;
	sub.data = malloc((sub.end-sub.start+1)*sizeof(double*));
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		sub.data[ii-sub.start] = malloc(ntheta*sizeof(double));
		for (ik=0; ik<ntheta; ik++)
		{
			sub.data[ii-sub.start][ik] = pol[ii][k][ik];
		}
	}

	bounds = malloc(3*sizeof(mp_par));

	for (ii=0; ii<3; ii++)
	{
		bounds[ii].fixed = 0;
		bounds[ii].side = 0;
		bounds[ii].deriv_debug = 0;
		bounds[ii].relstep = 0.0;
		bounds[ii].step = 0.0;
		bounds[ii].limited[0] = bounds[ii].limited[1] = 0;
	}
	bounds[0].limited[0] = 1;
	bounds[0].limits[0] = 0.0;
	bounds[2].limited[0] = 1;
	bounds[2].limits[0] = 0.0;

	mpconf->ftol = p->ftol;
	mpconf->xtol = p->xtol;
	mpconf->gtol = p->gtol;
	mpconf->covtol = 1e-10;
	mpconf->maxiter = p->niter;
	mpconf->nofinitecheck = 1;

	mpreturn = 0;
	mpreturn = mpfit(&funk_back, ntheta*(sub.end-sub.start+1), 3, 
					param, bounds, mpconf, &sub, mpres2);
	
	*amp = param[0];
	*cutoff = param[1];
	*power = param[2];

	free(param);
	free(bounds);
	free(mpres2);
	free(mpconf);
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		free(sub.data[ii-sub.start]);
	}
	free(sub.data);
	return 0;
}

int funk_back (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	int istw, iendw, num, iw, itht, ik;
	double back;

	sub = (struct kslice*) private;
	istw = sub->start;
	iendw = sub->end;

	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		back = p[0]/(1.+pow((iw*sub->delta_nu)/p[1],p[2]));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			deviates[num] = 0.0;
			if (sub->data[iw-istw][itht] > 0.0)
			{
				deviates[num] = (back - sub->data[iw-istw][itht])/back;
			}
			num++;
		}
	}
	return 0;
}
