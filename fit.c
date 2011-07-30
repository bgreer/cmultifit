#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "mpfit.h"
#include "header.h"

int fit_peak (float *freq, float *amp, float *width, float*** pol, float*** noise, float delta_nu, float delta_k, int nnu, int ntheta, int k)
{
	int ii, ij, mpreturn;
	double *param, *xerror;
	mp_par *bounds;
	mp_result *mpres2;
	mp_config *mpconf;
	struct kslice sub;

	param = (double*) malloc(4*sizeof(double));
	param[0] = *freq;
	param[1] = *amp*40.;
	param[2] = *width;
	param[3] = *amp*0.01;

	mpres2 = (mp_result*) calloc(1,sizeof(mp_result));
	mpconf = (mp_config*) calloc(1, sizeof(mp_config));

	/* Load subsection data */
	sub.start = (param[0] - 1.0*param[2])/delta_nu;
	sub.end = (param[0] + 1.0*param[2])/delta_nu;
	if (sub.start < 0) sub.start = 0;
	if (sub.end >= nnu) sub.end = nnu-1;
	sub.ntheta = ntheta;
	sub.delta_nu = delta_nu;
	sub.k = k*delta_k;
	sub.data = malloc((sub.end-sub.start+1)*sizeof(float*));
	sub.noise = malloc((sub.end-sub.start+1)*sizeof(float*));
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		sub.data[ii-sub.start] = malloc(ntheta*sizeof(float));
		memcpy(sub.data[ii-sub.start], pol[ii][k], ntheta*sizeof(float));
		sub.noise[ii-sub.start] = malloc(ntheta*sizeof(float));
		memcpy(sub.noise[ii-sub.start], noise[ii][k], ntheta*sizeof(float));
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
	bounds[2].limited[0] = bounds[2].limited[1] = 1;
	bounds[2].limits[0] = param[2]/3.;
	bounds[2].limits[1] = param[2]*3.;

	mpconf->ftol = 1e-8;
	mpconf->xtol = 1e-4;
	mpconf->gtol = 1e-8;
	mpconf->covtol = 1e-10;
	mpconf->maxiter = 200;

	mpres2->xerror = xerror;
	mpreturn = 1;
	mpreturn = mpfit(&funk_single, ntheta*(sub.end-sub.start+1), 4, 
					param, bounds, mpconf, &sub, mpres2);

	*freq = param[0];
	*amp = param[1]*0.4;/**2000./param[0];*/
	*width = param[2]*0.4;

	free(param);
	free(bounds);
	free(xerror);
	free(mpres2);
	free(mpconf);

	FILE *fp;
	float den;
	fp = fopen("debug2", "w");
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		den = ii*delta_nu - *freq;
		den = den*den + (*width)*(*width)/4.0;
		fprintf(fp, "%f\t%f\t%f\n", ii*delta_nu, param[3]+(*amp)*(*width)/(2.0*den), sub.data[ii-sub.start][0]);
	}
	fclose(fp);
	
	return 0;
}

int funk_single (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	int istw, iendw, num, iw, itht, ik;
	double akt, w1, twopi, tht, den, x2, err;

	sub = (struct kslice*) private;
	istw = sub->start;
	iendw = sub->end;
	akt = sub->k;
	twopi = 6.28318530717958647692528676655901;

	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = twopi*itht/sub->ntheta;
			deviates[num] = 0.0;
			if (sub->data[iw-istw][itht] > 0.0 && !isnan(sub->data[iw-istw][itht]))
			{
				den = w1 - p[0];
				den = den*den + p[2]*p[2]/4.0;
				x2 = p[1]*p[2]/(2.0*den) + p[3];
				deviates[num] = x2 - sub->data[iw-istw][itht];
				deviates[num] /= sub->noise[iw-istw][itht];
			}
			num++;
		}
	}
	return 0;
}

int fit_back (double* amp, double* cutoff, double* power, float*** pol, float*** noise, float delta_nu, int nnu, int ntheta, int k)
{
	int ii, ij, mpreturn;
	double *param;
	mp_par *bounds;
	mp_result *mpres2;
	mp_config *mpconf;
	struct kslice sub;

	param = (double*) malloc(3*sizeof(double));
	param[0] = *amp;
	param[1] = *cutoff;
	param[2] = *power;

	mpres2 = (mp_result*) calloc(1,sizeof(mp_result));
	mpconf = (mp_config*) calloc(1, sizeof(mp_config));

	/* Load subsection data */
	sub.start = 0;
	sub.end = (1000.)/delta_nu;
	if (sub.end >= nnu) sub.end = nnu-1;
	sub.ntheta = ntheta;
	sub.delta_nu = delta_nu;
	sub.k = 0.0;
	sub.data = malloc((sub.end-sub.start+1)*sizeof(float*));
	sub.noise = malloc((sub.end-sub.start+1)*sizeof(float*));
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		sub.data[ii-sub.start] = malloc(ntheta*sizeof(float));
		memcpy(sub.data[ii-sub.start], pol[ii][k], ntheta*sizeof(float));
		sub.noise[ii-sub.start] = malloc(ntheta*sizeof(float));
		memcpy(sub.noise[ii-sub.start], noise[ii][k], ntheta*sizeof(float));
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

	mpconf->ftol = 1e-8;
	mpconf->xtol = 1e-8;
	mpconf->gtol = 1e-8;
	mpconf->covtol = 1e-10;
	mpconf->maxiter = 200;

	mpreturn = 1;
	mpreturn = mpfit(&funk_back, ntheta*(sub.end-sub.start+1), 3, 
					param, bounds, mpconf, &sub, mpres2);

	*amp = param[0];
	*cutoff = param[1];
	*power = param[2];
	
	FILE *fp;
	float den;
	fp = fopen("debug2", "w");
	for (ii=sub.start; ii<=sub.end; ii++)
	{
		fprintf(fp, "%f\t%f\t%f\n", ii*delta_nu, *amp/(1.+pow(*cutoff*ii*delta_nu,*power)), sub.data[ii-sub.start][0]);
	}
	fclose(fp);
	


	free(param);
	free(bounds);
	free(mpres2);
	free(mpconf);
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
		back = p[0]/(1.+pow(p[1]*iw*sub->delta_nu,p[2]));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			deviates[num] = (back - sub->data[iw-istw][itht])/sub->noise[iw-istw][itht];
			num++;
		}
	}
	return 0;
}
