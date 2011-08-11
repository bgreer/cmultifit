#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

void usage (char* name)
{
	printf("Usage: \n%s paramfile\n", name);
}

int main (int argc, char* argv[])
{
	FILE *fpmodel, *fpout, *fpback;
	struct params par;
	int status = 0;
	int ii, ij, ik;
	int ntheta, nk, nnu;
	double delta_nu, delta_k;
	double ***pol, ***noise, *norm;
	double **freq, **amp, **width;
	float readfreq, readamp, readwidth, readk;
	int *numridges;

	int mpreturn, index;
	struct kslice subsection;
	double* param;
	mp_result *mpres;
	mp_config *mpconf;
	mp_par *bounds;
	double *xerror, *covar;
	mpres = malloc(sizeof(mp_result));
	mpconf = malloc(sizeof(mp_config));
	double sum, sum2, sum3;

	if (argc < 2)
	{
		usage(argv[0]);
		return EXIT_FAILURE;
	}

	/* Read param file and set run parameters */
	read_param_file(argv[1], &par);

	/* Open FITS file and read keys */
	if (read_fits_file(&pol, &noise, &par, &ntheta, &nk, &nnu, &delta_k, &delta_nu)==EXIT_FAILURE)
		return EXIT_FAILURE;
	thtarr = malloc(ntheta*sizeof(double));
	thtpow = malloc(ntheta*sizeof(double));

	/* Normalize spectrum */
	norm = (double*) calloc(nk, sizeof(double));
	/*normalize(&pol, &norm, nnu, nk, ntheta);*/

	/* TODO: Check kstart and kend */
	if (par.kend >= nk) par.kend = nk-1;

/* Read model file */
	if (!par.silent) printf("Reading model from %s\n", par.modelfname);
	numridges = (int*) malloc(nk*sizeof(int));
	freq = (double**) malloc(nk*sizeof(double*));
	amp = (double**) malloc(nk*sizeof(double*));
	width = (double**) malloc(nk*sizeof(double*));
	for (ii=0; ii<nk; ii++)
		numridges[ii] = 0;
	/* Run through once to count ridges */
	fpmodel = fopen(par.modelfname, "r");
	if (fpmodel==NULL)
	{
		printf("ERROR: could not open model file: %s\n", par.modelfname);
		return EXIT_FAILURE;
	}
	while (fscanf(fpmodel, "%d\t%e\f%e\t%e\t%e\n", &ii, &readk, &readfreq, &readamp, &readwidth) != EOF)
	{
		if (ii>=0 && ii<nk)
		{
			numridges[ii]++;
		} else {
			printf("ERROR: k-value out of range in model: %d\n", ii);
			return EXIT_FAILURE;
		}
	}
	rewind(fpmodel);
	/* allocate enough for all ridges, then load them */
	for (ii=0; ii<nk; ii++)
	{
		if (numridges[ii]>0)
		{
			freq[ii] = (double*) malloc(numridges[ii]*sizeof(double));
			amp[ii] = (double*) malloc(numridges[ii]*sizeof(double));
			width[ii] = (double*) malloc(numridges[ii]*sizeof(double));
		}
		for (ij=0; ij<numridges[ii]; ij++)
		{
			fscanf(fpmodel, "%d\t%e\t%e\t%e\t%e\n", &ik, &readk, &readfreq, &readamp, &readwidth);
			if (readfreq<0.0 || readamp<0.0 || readwidth<0.0)
			{
				printf("ERROR: Invalid parameters in model, k=%d\n", ii);
				return EXIT_FAILURE;
			}
			freq[ii][ij] = readfreq;
			amp[ii][ij] = readamp;
			width[ii][ij] = readwidth;
		}
	}
	fclose(fpmodel);

/* Pre-compute noise cube */
	if (par.chiweight == WEIGHT_NOISE)
	{
		if (!par.silent) printf("Pre-computing noise.\n");
		switch (par.noisemode)
		{
			case NOISE_CONST:
				compute_noise_const(pol, noise, nnu, nk, ntheta);
				break;
			case NOISE_SMOOTH:
				compute_noise_smooth(pol, noise, nnu, nk, ntheta, (800.)/delta_nu);
				break;
			case NOISE_WAVELET:
				compute_noise_wavelet(pol, noise, nnu, nk, ntheta);
				break;
		}
	}

/* Open output file */
	fpout = fopen(par.outfname, "w");
	if (fpout==NULL)
	{
		printf("ERROR: could not open output file: %s\n", par.outfname);
		return EXIT_FAILURE;
	}

/* Open background output file if needed */
	if (par.backfname)
	{
		fpback = fopen(par.backfname, "w");
		if (fpback==NULL)
		{
			printf("ERROR: could not open output file: %s\n", par.backfname);
			return EXIT_FAILURE;
		}
	}

	if (!par.silent) printf("\nBeginning optimization, output file '%s'\n", par.outfname);

/* Enter loop over all k */
	for (ij=par.kstart; ij<=par.kend; ij++)
	{
		if (numridges[ij]>0)
		{
			if (!par.silent) printf("Starting fit of %d ridges at k=%d\n", numridges[ij], ij);
			param = (double*) malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			bounds = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(mp_par));
			xerror = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			if (par.covarfname)
				covar = malloc(((numridges[ij]*NPEAK+NBACK)*(numridges[ij]*NPEAK+NBACK))*sizeof(double));

			/* Do rough fit of single peaks */
			/* TODO: fix this */
			if (!par.silent) printf("\tDoing single ridge estimates.\n");
			for (ii=0; ii<numridges[ij]; ii++)
			{
				fit_peak(&par, &(freq[ij][ii]), &(amp[ij][ii]), &(width[ij][ii]), pol, noise, 
					delta_nu, delta_k, nnu, ntheta, ij);
			}
			
			/* Do rough fit of background */
			if (!par.silent) printf("\tFitting background at low frequency.\n");
			fit_back(&par, &(param[numridges[ij]*NPEAK]), &(param[numridges[ij]*NPEAK+1]), 
				&(param[numridges[ij]*NPEAK+2]), pol, noise, delta_nu, nnu, ntheta, ij);
			
			/* Set multifit contraints */
			for (ii=0; ii<numridges[ij]; ii++)
			{
				for (ik=0; ik<NPEAK; ik++)
				{
					bounds[ii*NPEAK+ik].fixed = 0;
					bounds[ii*NPEAK+ik].side = 0;
					bounds[ii*NPEAK+ik].deriv_debug = 0;
					bounds[ii*NPEAK+ik].relstep = 0.0;
					bounds[ii*NPEAK+ik].step = 0.0;
				}

				/* set frequency */
				param[ii*NPEAK] = freq[ij][ii];
				bounds[ii*NPEAK].limited[0] = bounds[ii*NPEAK].limited[1] = 1;
				if (ii==0)
					bounds[ii*NPEAK].limits[0] = 0.0;
				else 
					bounds[ii*NPEAK].limits[0] = 0.5*(freq[ij][ii] + freq[ij][ii-1]);
				if (ii==numridges[ij]-1)
					bounds[ii*NPEAK].limits[1] = nnu*delta_nu;
				else
					bounds[ii*NPEAK].limits[1] = 0.5*(freq[ij][ii] + freq[ij][ii+1]);

				/* set amplitude */
				param[ii*NPEAK+1] = amp[ij][ii];
				bounds[ii*NPEAK+1].limited[0] = 1;
				bounds[ii*NPEAK+1].limited[1] = 0;
				bounds[ii*NPEAK+1].limits[0] = 0.0; /* dont go below 0 */

				/* set width */
				param[ii*NPEAK+2] = width[ij][ii];
				bounds[ii*NPEAK+2].limited[0] = bounds[ii*NPEAK+2].limited[1] = 1;
				bounds[ii*NPEAK+2].limits[0] = width[ij][ii]/4.0;
				bounds[ii*NPEAK+2].limits[1] = width[ij][ii]*4.0;

				/* set velocities */
				param[ii*NPEAK+3] = 1.0;
				bounds[ii*NPEAK+3].limited[0] = bounds[ii*NPEAK+3].limited[1] = 1;
				bounds[ii*NPEAK+3].limits[0] = -1000.0;
				bounds[ii*NPEAK+3].limits[1] = 1000.0;
				param[ii*NPEAK+4] = 1.0;
				bounds[ii*NPEAK+4].limited[0] = bounds[ii*NPEAK+4].limited[1] = 1;
				bounds[ii*NPEAK+4].limits[0] = -1000.0;
				bounds[ii*NPEAK+4].limits[1] = 1000.0;

				/* set theta variation params */
				param[ii*NPEAK+5] = 0.1;
				bounds[ii*NPEAK+5].limited[0] = bounds[ii*NPEAK+5].limited[1] = 1;
				bounds[ii*NPEAK+5].limits[0] = -1.0;
				bounds[ii*NPEAK+5].limits[1] = 1.0;
				param[ii*NPEAK+6] = PI/2.0;
				bounds[ii*NPEAK+6].limited[0] = bounds[ii*NPEAK+6].limited[1] = 0;
			}
			for (ii=0; ii<NBACK; ii++)
			{
				bounds[numridges[ij]*NPEAK+ii].fixed = 0;
				bounds[numridges[ij]*NPEAK+ii].side = 0;
				bounds[numridges[ij]*NPEAK+ii].deriv_debug = 0;
				bounds[numridges[ij]*NPEAK+ii].relstep = 0.0;
				bounds[numridges[ij]*NPEAK+ii].step = 0.0;
			}
			bounds[numridges[ij]*NPEAK].limited[0] = 1;
			bounds[numridges[ij]*NPEAK].limited[1] = 0;
			bounds[numridges[ij]*NPEAK].limits[0] = 0.0;
			bounds[numridges[ij]*NPEAK+1].limited[0] = bounds[numridges[ij]*NPEAK+1].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+2].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+2].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+2].limits[0] = 0.0;
			param[numridges[ij]*NPEAK+3] = 0.1;
			bounds[numridges[ij]*NPEAK+3].limited[0] = bounds[numridges[ij]*NPEAK+3].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+3].limits[0] = -1.0;
			bounds[numridges[ij]*NPEAK+3].limits[1] = 1.0;
			param[numridges[ij]*NPEAK+4] = PI/2.0;
			bounds[numridges[ij]*NPEAK+4].limited[0] = bounds[numridges[ij]*NPEAK+4].limited[1] = 0;

			param[numridges[ij]*NPEAK+5] = 0.1*param[numridges[ij]*NPEAK];
			bounds[numridges[ij]*NPEAK+5].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+5].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+5].limits[0] = 0.0;
			param[numridges[ij]*NPEAK+6] = 1500.;
			bounds[numridges[ij]*NPEAK+6].limited[0] = bounds[numridges[ij]*NPEAK+6].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+6].limits[0] = 200.;
			bounds[numridges[ij]*NPEAK+6].limits[1] = nnu*delta_nu;
			param[numridges[ij]*NPEAK+7] = 1000.;
			bounds[numridges[ij]*NPEAK+7].limited[0] = bounds[numridges[ij]*NPEAK+7].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+7].limits[0] = 500.0;
			bounds[numridges[ij]*NPEAK+7].limits[1] = 10000.;

			subsection.par = &par;
			
			/* Load klice struct for passing to function */
			subsection.start = 0;
			if (subsection.start < 0) subsection.start = 0;
			
			subsection.end = (param[(numridges[ij]-1)*NPEAK]+1.*param[(numridges[ij]-1)*NPEAK+2])
								/delta_nu;
			if (subsection.end > nnu-1) subsection.end = nnu-1;

			subsection.delta_nu = delta_nu;
			subsection.ntheta = ntheta;
			subsection.k = (ij+1)*delta_k;
			subsection.data = malloc((subsection.end-subsection.start+1)*sizeof(double*));
			subsection.noise = malloc((subsection.end-subsection.start+1)*sizeof(double*));
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				subsection.data[ii-subsection.start] = malloc(ntheta*sizeof(double));
				subsection.noise[ii-subsection.start] = malloc(ntheta*sizeof(double));
				for (ik=0; ik<ntheta; ik++)
				{
					subsection.data[ii-subsection.start][ik] = pol[ii][ij][ik];
					subsection.noise[ii-subsection.start][ik] = noise[ii][ij][ik];
				}
			}

			/* Set optimization parameters */
			mpconf->ftol = par.ftol;
			mpconf->xtol = par.xtol;
			mpconf->gtol = par.gtol;
			mpconf->covtol = 1e-14;
			mpconf->maxiter = par.niter;
			mpconf->nofinitecheck = 1;

			mpres->xerror = xerror;
			if (par.covarfname) mpres->covar = covar;
			mpreturn = 1;
/* Perform optimization */
			if (!par.silent) printf("\tDoing multifit.\n");

			mpreturn = mpfit(&funk, 
					ntheta*(subsection.end-subsection.start+1), 
					numridges[ij]*NPEAK+NBACK, 
					param, 
					bounds, 
					mpconf, 
					&subsection, 
					mpres);
		
			/* Check return value of mpfit */
			switch (mpreturn)
			{
				case MP_OK_CHI:
					printf("\tMPFIT convergence in chi-square\n");break;
				case MP_OK_PAR:
					printf("\tMPFIT convergence in parameter value\n");break;
				case MP_OK_DIR:
					printf("\tMPFIT convergence in orthogonality\n");break;
				case MP_MAXITER:
					printf("\tMPFIT max iterations reached\n");break;
				case MP_FTOL:
					printf("\tMPFIT ftol criteria reached\n");break;
				case MP_XTOL:
					printf("\tMPFIT xtol criteria reached\n");break;
				case MP_GTOL:
					printf("\tMPFIT gtol criteria reached\n");break;

				case MP_ERR_INPUT:
					printf("\tMPFIT ERROR: input parameter error\n");break;
				case MP_ERR_NAN:
					printf("\tMPFIT ERROR: function returned nan\n");break;
				case MP_ERR_MEMORY:
					printf("\tMPFIT ERROR: memory allocation error\n");break;
				case MP_ERR_INITBOUNDS:
					printf("\tMPFIT ERROR: guesses inconsistent with bounds\n");break;
				case MP_ERR_BOUNDS:
					printf("\tMPFIT ERROR: inconsistent bounds\n");break;
				case MP_ERR_PARAM:
					printf("\tMPFIT ERROR: input parameter error\n");break;
				case MP_ERR_DOF:
					printf("\tMPFIT ERROR: too few degrees of freedom\n");break;

				default:
					printf("\tMPFIT return = %d\n", mpreturn);break;
			}

			/* Output optimization details */
			if (!par.silent) 
			{
				printf("\tStarting chi2 = %e\n", 
					mpres->orignorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tFinal chi2 = %e\n", 
					mpres->bestnorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tNumber of iterations = %d\n", mpres->niter);
				printf("\tNumber of function evals = %d\n", mpres->nfev);
				printf("\tNumber of pegged parameters = %d\n", mpres->npegged);
			}

			/* Post-process fit parameters */
			for (ii=0; ii<numridges[ij]; ii++)
			{
				if (param[ii*NPEAK+5] < 0.0)
				{
					param[ii*NPEAK+6] += PI/2.0;
					param[ii*NPEAK+5] = fabs(param[ii*NPEAK+5]);
				}
				param[ii*NPEAK+6] = fmod(param[ii*NPEAK+6], PI);
				if (param[ii*NPEAK+6] < 0.0)
					param[ii*NPEAK+6] += PI;
			}
			if (param[numridges[ij]*NPEAK+3] < 0.0)
			{
				param[numridges[ij]*NPEAK+4] += PI/2.0;
				param[numridges[ij]*NPEAK+3] = fabs(param[numridges[ij]*NPEAK+3]);
			}
			param[numridges[ij]*NPEAK+4] = fmod(param[numridges[ij]*NPEAK+4], PI);
			if (param[numridges[ij]*NPEAK+4] < 0.0)
				param[numridges[ij]*NPEAK+4] += PI;

			/* Print fit debug */
			if (par.debugfname) output_debug(&par, pol, noise, ntheta, nk, nnu, ij, ntheta*(subsection.end-subsection.start+1), 
					numridges[ij]*NPEAK+NBACK, 
					param, delta_nu, delta_k);
		
			/* Print covariance matrix */
			if (par.covarfname) output_covar(covar, numridges[ij], &par);

			/* Output fit to file if valid */
			if (mpreturn!=MP_MAXITER && mpreturn > 0 && mpreturn!=3)
			{
				for (ii=0; ii<numridges[ij]; ii++)
				{
					if (param[ii*NPEAK+2] < 7.5)
					{
						printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
						printf("Width too small (%f)\n", param[ii*NPEAK+2]);
					} else if (param[ii*NPEAK+2] > 800.) {
						printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
						printf("Width too large (%f)\n", param[ii*NPEAK+2]);
					} else if (abs(param[ii*NPEAK]-freq[ij][ii])/freq[ij][ii] >= 0.2) {
						printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
						printf("Frequency drifted too much from model\n");
					} else if (xerror[ii*NPEAK]<=0.0 || xerror[ii*NPEAK+1]<=0.0 || 
							xerror[ii*NPEAK+2]<=0.0 || xerror[ii*NPEAK+3]<=0.0 || 
							xerror[ii*NPEAK+4]<=0.0) {
						printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
						printf("Zero error\n");
					} else {
						fprintf(fpout, "%d\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\n", 
							ij, 
							ij*delta_k,
							param[ii*NPEAK]/(ij*delta_k), 
							param[ii*NPEAK], xerror[ii*NPEAK+3], 
							param[ii*NPEAK+1], xerror[ii*NPEAK+1], 
							param[ii*NPEAK+2], xerror[ii*NPEAK+2], 
							param[ii*NPEAK+3], xerror[ii*NPEAK+3], 
							param[ii*NPEAK+4], xerror[ii*NPEAK+4],
							param[ii*NPEAK+5], xerror[ii*NPEAK+5], 
							param[ii*NPEAK+6], xerror[ii*NPEAK+6]
							);
					}
				}
				fflush(fpout);
			
				/* Output background params to separate file */
				if (par.backfname)
				{
					fprintf(fpback, "%d\t%f", ij, ij*delta_k);
					for (ik=0; ik<NBACK; ik++)
						fprintf(fpback, "\t%e", param[numridges[ij]*NPEAK+ik]);
					fprintf(fpback, "\n");
					fflush(fpback);
				}

			} else {
				if (!par.silent) printf("\tFit deemed invalid, nothing printed to file\n");
			}

			/* Free some memory for next k */
			free(param);
			free(bounds);
			free(xerror);
			if (par.covarfname) free(covar);
		}
	}
	fclose(fpout);
	if (par.backfname) fclose(fpback);
	return EXIT_SUCCESS;
}

double model (int numridges, double nu, double k, double theta, double* p)
{
	int ii;
	double ret, den;

	ret = 0.0;
	for (ii=0; ii<numridges; ii++)
	{
		den = nu + k*(p[ii*NPEAK+3]*cos(theta)+p[ii*NPEAK+4]*sin(theta))/6.28318531 - p[ii*NPEAK];
		den = den*den + p[ii*NPEAK+2]*p[ii*NPEAK+2]/4.0;
		ret += p[ii*NPEAK+1]*p[ii*NPEAK+2]/(2.0*den);
	}
	ret += p[numridges*NPEAK]/(1.0 + pow(nu*p[numridges*NPEAK+1], p[numridges*NPEAK+2]));
	ret += p[numridges*NPEAK+3]*p[numridges*NPEAK+5]/(pow(nu-p[numridges*NPEAK+4],2.) + pow(0.5*p[numridges*NPEAK+5],2.));
/*	ret += p[numridges*5+3]*p[numridges*5+5]/((nu-p[numridges*5+4])*(nu-p[numridges*5+4])+p[numridges*5+5]*p[numridges*5+5]/4.);
*/	return ret;
}

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

/* Normalize spectrum by average at constant k */
void normalize (double ****spec, double** norm, int nnu, int nk, int ntheta)
{
	int ii, ij, ik;
	double sum;

	for (ij=0; ij<nk; ij++)
	{
		sum = 0.0;
		for (ii=0; ii<nnu; ii++)
			for (ik=0; ik<ntheta; ik++)
				sum += (*spec)[ii][ij][ik];
		(*norm)[ij] = sum/(nnu*ntheta);
		for (ii=0; ii<nnu; ii++)
			for (ik=0; ik<ntheta; ik++)
				(*spec)[ii][ij][ik] /= (*norm)[ij];
	}
}
