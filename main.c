#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"

void usage (char* name)
{
	printf("Usage:\n%s spectrum model output\n", name);
}

int main (int argc, char* argv[])
{
	fitsfile *fptr;
	FILE *fpmodel, *fpout, *fpdebug;
	int status = 0;
	int silent, ii, ij, ik;
	int ntheta, nk, nnu;
	float delta_nu, delta_k;
	float ***pol, *buff, ***noise;
	long coords[3];
	float **freq, **amp, **width;
	float readfreq, readamp, readwidth, readk;
	int *numridges;

	int mpreturn, index, numback;
	struct kslice subsection;
	double* param;
	mp_result *mpres;
	mp_config *mpconf;
	mp_par *bounds;
	double *xerror, *covar;
	mpres = malloc(sizeof(mp_result));
	mpconf = malloc(sizeof(mp_config));
	double sum, sum2, sum3;

	if (argc < 4)
	{
		usage(argv[0]);
		return EXIT_FAILURE;
	}
	silent = 0;
	numback = 6;

/* Open FITS file and read keys */
	if (!silent) printf("Reading FITS file: %s\n", argv[1]);
	fits_open_file(&fptr, argv[1], READONLY, &status);
	if (status)
	{
		printf("Error opening FITS file: %s\n", argv[1]);
		fits_report_error(stdout, status);
		return EXIT_FAILURE;
	}
	fits_read_key(fptr, TINT, "NAXIS1", &ntheta, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", &nk, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS3", &nnu, NULL, &status);
	if (status)
	{
		printf("Could not find FITS dimensions.\n");
		printf("Make sure the following keys are set: NAXIS1, NAXIS2, NAXIS3\n");
		return EXIT_FAILURE;
	}
	if (!silent) printf("FITS dimensions:\n\ttheta:\t%d\n\tk:\t%d\n\tnu:\t%d\n", ntheta, nk, nnu);
	fits_read_key(fptr, TFLOAT, "DELTA_K", &delta_k, NULL, &status);
	fits_read_key(fptr, TFLOAT, "DELTA_NU", &delta_nu, NULL, &status);
	if (status)
	{
		printf("Could not find FITS keys DELTA_NU and/or DELTA_K.\n");
		return EXIT_FAILURE;
	}
	if (!silent) printf("\tDELTA_NU:\t%f\n\tDELTA_K:\t%f\n", delta_nu, delta_k);

/* Allocate memory for entire FITS data cube */
	if (!silent) printf("Allocating %ld bytes for data cube.\n", sizeof(float)*ntheta*nk*nnu);
	if (!silent) printf("\t With %ld bytes of overhead.\n", sizeof(float*)*nk*nnu + sizeof(float**)*nnu);

	pol = (float***) malloc(nnu*sizeof(float**));
	noise = (float***) malloc(nnu*sizeof(float**));
	for (ii=0; ii<nnu; ii++)
	{
		pol[ii] = (float**) malloc(nk*sizeof(float*));
		noise[ii] = (float**) malloc(nk*sizeof(float*));
		for (ij=0; ij<nk; ij++)
		{
			pol[ii][ij] = (float*) malloc(ntheta*sizeof(float));
			noise[ii][ij] = (float*) malloc(ntheta*sizeof(float));
			if (pol[ii][ij]==NULL)
			{
				printf("Error allocating memory.\n");
				return EXIT_FAILURE;
			}
		}
	}
	buff = (float*) malloc(ntheta*sizeof(float));

/* Read FITS file into array */
	coords[0] = 1L;
	if (!silent) printf("Reading data cube into memory.\n");
	for (coords[2]=1; coords[2]<=nnu; coords[2]++)
	{
		for (coords[1]=1; coords[1]<=nk; coords[1]++)
		{
			fits_read_pix(fptr, TFLOAT, coords, ntheta, NULL, buff, NULL, &status);
			for (ik=0; ik<ntheta; ik++)
			{
				if (buff[ik] < 0.0 || isnan(buff[ik])) buff[ik] = 0.0;
			}
			memcpy(pol[coords[2]-1][coords[1]-1], buff, ntheta*sizeof(float));
		}
	}
	fits_close_file(fptr, &status);

/* Read model file */
	if (!silent) printf("Reading model from %s\n", argv[2]);
	numridges = (int*) malloc(nk*sizeof(int));
	freq = (float**) malloc(nk*sizeof(float*));
	amp = (float**) malloc(nk*sizeof(float*));
	width = (float**) malloc(nk*sizeof(float*));
	for (ii=0; ii<nk; ii++)
		numridges[ii] = 0;
	/* Run through once to count ridges */
	fpmodel = fopen(argv[2], "r");
	if (fpmodel==NULL)
	{
		printf("ERROR: could not open model file: %s\n", argv[2]);
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
			freq[ii] = (float*) malloc(numridges[ii]*sizeof(float));
			amp[ii] = (float*) malloc(numridges[ii]*sizeof(float));
			width[ii] = (float*) malloc(numridges[ii]*sizeof(float));
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
	if (!silent) printf("Pre-computing noise.\n");
	compute_noise_wavelet(pol, noise, nnu, nk, ntheta, (800.)/delta_nu);

/* Open output file */
	fpout = fopen(argv[3], "w");
	if (fpout==NULL)
	{
		printf("ERROR: could not open output file: %s\n", argv[3]);
		return EXIT_FAILURE;
	}

	printf("Beginning optimization.\nOutput file '%s' will be in the format:\n", argv[3]);
	printf("k-bin  k  w/k  freq  err  amp  err  width  err  ux  err uy  err\n");

/* Enter loop over all k */
/*	for (ij=0; ij<nk; ij++)*/
	ij = 14;
	{
		if (numridges[ij]>0)
		{
			if (!silent) printf("Starting fit of %d ridges at k=%d\n", numridges[ij], ij);
			param = (double*) malloc((numridges[ij]*7+numback)*sizeof(double));
			bounds = malloc((numridges[ij]*7+numback)*sizeof(mp_par));
			xerror = malloc((numridges[ij]*7+numback)*sizeof(double));
			covar = malloc(((numridges[ij]*7+numback)*(numridges[ij]*7+numback))*sizeof(double));
			/* Do rough fit of single peaks */
			if (!silent) printf("\tDoing single ridge estimates.\n");
			for (ii=0; ii<numridges[ij]; ii++)
			{
				fit_peak(&(freq[ij][ii]), &(amp[ij][ii]), &(width[ij][ii]), pol, noise, 
					delta_nu, delta_k, nnu, ntheta, ij);
			}

			for (ii=0; ii<numridges[ij]; ii++)
			{
				for (ik=0; ik<7; ik++)
				{
					bounds[ii*7+ik].fixed = 0;
					bounds[ii*7+ik].side = 0;
					bounds[ii*7+ik].deriv_debug = 0;
					bounds[ii*7+ik].relstep = 0.0;
					bounds[ii*7+ik].step = 0.0;
				}

				/* set frequency */
				param[ii*7] = freq[ij][ii];
				bounds[ii*7].limited[0] = bounds[ii*7].limited[1] = 1;
				if (ii==0)
					bounds[ii*7].limits[0] = 0.0;
				else 
					bounds[ii*7].limits[0] = 0.5*(freq[ij][ii] + freq[ij][ii-1]);
				if (ii==numridges[ij]-1)
					bounds[ii*7].limits[1] = nnu*delta_nu;
				else
					bounds[ii*7].limits[1] = 0.5*(freq[ij][ii] + freq[ij][ii+1]);

				/* set amplitude */
				param[ii*7+1] = amp[ij][ii];
				bounds[ii*7+1].limited[0] = 1;
				bounds[ii*7+1].limited[1] = 0;
				bounds[ii*7+1].limits[0] = 0.0; /* dont go below 0 */

				/* set width */
				param[ii*7+2] = width[ij][ii];
				bounds[ii*7+2].limited[0] = bounds[ii*7+2].limited[1] = 1;
				bounds[ii*7+2].limits[0] = width[ij][ii]/4.0;
				bounds[ii*7+2].limits[1] = width[ij][ii]*4.0;

				/* set velocities */
				param[ii*7+3] = 1.0;
				bounds[ii*7+3].limited[0] = bounds[ii*7+3].limited[1] = 1;
				bounds[ii*7+3].limits[0] = -1000.0;
				bounds[ii*7+3].limits[1] = 1000.0;
				param[ii*7+4] = 1.0;
				bounds[ii*7+4].limited[0] = bounds[ii*7+4].limited[1] = 1;
				bounds[ii*7+4].limits[0] = -1000.0;
				bounds[ii*7+4].limits[1] = 1000.0;

				/* set theta variation params */
				param[ii*7+5] = 0.0;
				bounds[ii*7+5].limited[0] = 0;
				bounds[ii*7+5].limited[1] = 0;
				param[ii*7+6] = 0.0;
				bounds[ii*7+6].limited[0] = bounds[ii*7+6].limited[1] = 0;
			}
			for (ii=0; ii<numback; ii++)
			{
				bounds[numridges[ij]*7+ii].fixed = 0;
				bounds[numridges[ij]*7+ii].side = 0;
				bounds[numridges[ij]*7+ii].deriv_debug = 0;
				bounds[numridges[ij]*7+ii].relstep = 0.0;
				bounds[numridges[ij]*7+ii].step = 0.0;
			}
			param[numridges[ij]*7] = 0.5;
			bounds[numridges[ij]*7].limited[0] = 1;
			bounds[numridges[ij]*7].limited[1] = 0;
			bounds[numridges[ij]*7].limits[0] = 0.0;
			param[numridges[ij]*7+1] = 0.0033;
			bounds[numridges[ij]*7+1].limited[0] = bounds[numridges[ij]*7+1].limited[1] = 0;
			param[numridges[ij]*7+2] = 2.5;
			bounds[numridges[ij]*7+2].limited[0] = 1;
			bounds[numridges[ij]*7+2].limited[1] = 0;
			bounds[numridges[ij]*7+2].limits[0] = 0.0;

			param[numridges[ij]*7+3] = 2.;	
			bounds[numridges[ij]*7+3].limited[0] = 1;
			bounds[numridges[ij]*7+3].limited[1] = 0;
			bounds[numridges[ij]*7+3].limits[0] = 0.0;
			param[numridges[ij]*7+4] = 2000.;
			bounds[numridges[ij]*7+4].limited[0] = bounds[numridges[ij]*7+4].limited[1] = 1;
			bounds[numridges[ij]*7+4].limits[0] = 200.;
			bounds[numridges[ij]*7+4].limits[1] = nnu*delta_nu;
			param[numridges[ij]*7+5] = 1000.;
			bounds[numridges[ij]*7+5].limited[0] = bounds[numridges[ij]*7+5].limited[1] = 1;
			bounds[numridges[ij]*7+5].limits[0] = 500.0;
			bounds[numridges[ij]*7+5].limits[1] = 10000.;

			if (!silent) printf("\tFitting background at low frequency.\n");
			fit_back(&(param[numridges[ij]*7]), &(param[numridges[ij]*7+1]), 
				&(param[numridges[ij]*7+2]), pol, noise, delta_nu, nnu, ntheta, ij);
			/* Load klice struct for passing to function */
			subsection.start = 0;
			if (subsection.start < 0) subsection.start = 0;
			subsection.end = (param[(numridges[ij]-1)*7]+1.*param[(numridges[ij]-1)*7+2])/delta_nu;
			if (subsection.end > nnu-1) subsection.end = nnu-1;
			subsection.delta_nu = delta_nu;
			subsection.ntheta = ntheta;
			subsection.k = (ij+1)*delta_k;
			subsection.data = malloc((subsection.end-subsection.start+1)*sizeof(float*));
			subsection.noise = malloc((subsection.end-subsection.start+1)*sizeof(float*));
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				subsection.data[ii-subsection.start] = malloc(ntheta*sizeof(float));
				memcpy(subsection.data[ii-subsection.start], pol[ii][ij], ntheta*sizeof(float));
				subsection.noise[ii-subsection.start] = malloc(ntheta*sizeof(float));
				memcpy(subsection.noise[ii-subsection.start], noise[ii][ij], ntheta*sizeof(float));
			}

			/* Set optimization parameters */
			mpconf->ftol = 1e-8;
			mpconf->xtol = 1e-4;
			mpconf->gtol = 1e-8;
			mpconf->covtol = 1e-10;
			mpconf->maxiter = 200;
			mpconf->nofinitecheck = 1;

			mpres->xerror = xerror;
			mpres->covar = covar;
			mpreturn = 1;
/* Perform optimization */
			if (!silent) printf("\tDoing multifit.\n");
			mpreturn = mpfit(&funk, ntheta*(subsection.end-subsection.start+1), 
				numridges[ij]*7+numback, param, bounds, mpconf, &subsection, mpres);
		
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
			if (!silent) 
			{
				printf("\tStarting chi2 = %e\n", mpres->orignorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tFinal chi2 = %e\n", mpres->bestnorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tNumber of iterations = %d\n", mpres->niter);
				printf("\tNumber of function evals = %d\n", mpres->nfev);
				printf("\tNumber of pegged parameters = %d\n", mpres->npegged);
			}

			for (ii=0; ii<numback; ii++)
				printf("%f\n", param[numridges[ij]*7+ii]);

			fpdebug = fopen("debug", "w");
			for (ii=0; ii<nnu; ii++)
			{
				sum = sum2 = sum3 = 0.0;
				for (ik=0; ik<1; ik++)
				{
					sum = model(numridges[ij], ii*delta_nu, ij*delta_k, 6.28318531*(ik+1)/ntheta, param)/1.0;
					sum2 = pol[ii][ij][ik]/1.0;
					sum3 = noise[ii][ij][ik]/1.0;
				
				fprintf(fpdebug, "%f\t%e\t%e\t%e\t%e\t%e\n", ii*delta_nu, sum2, sum, 
					param[numridges[ij]*7]/(1.+pow(ii*delta_nu*param[numridges[ij]*7+1],
						param[numridges[ij]*7+2])), 
					param[numridges[ij]*7+3]*param[numridges[ij]*7+5]/
						(pow(ii*delta_nu-param[numridges[ij]*7+4], 2.) + 
							pow(0.5*param[numridges[ij]*7+5],2.)), 
					sum3);
				}
				fprintf(fpdebug, "");
			}
			fclose(fpdebug);

			for (ii=0; ii<numridges[ij]*7+numback; ii++)
			{
				for (ik=0; ik<numridges[ij]*7+numback; ik++)
				{
					printf("%d\t%d\t%e\n", ii, ik, log10(fabs(covar[ii*(numridges[ij]*7+numback)+ik]+1e-7)));
				}
				printf("\n");
			}

			/* Output fit to file if valid */
			if (mpreturn!=MP_MAXITER && mpreturn > 0 && mpreturn!=3)
			{
				for (ii=0; ii<numridges[ij]; ii++)
				{
					if (param[ii*7+2] > 7.5 && param[ii*7+2] < 665. && 
							abs(param[ii*7]-freq[ij][ii])/freq[ij][ii] < 0.2 && 
							param[ii*7+3] < 900. && param[ii*7+4] < 900. &&
							xerror[ii*7]>0.0 && xerror[ii*7+1]>0.0 && xerror[ii*7+2]>0.0 && 
							xerror[ii*7+3]>0.0 && xerror[ii*7+4]>0.0)
					{
						fprintf(fpout, "%d\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\n", 
							ij, 
							ij*delta_k,
							param[ii*7]/(ij*delta_k), 
							param[ii*7], xerror[ii*7+3], 
							param[ii*7+1], xerror[ii*7+1], 
							param[ii*7+2], xerror[ii*7+2], 
							param[ii*7+3], xerror[ii*7+3], 
							param[ii*7+4], xerror[ii*7+4],
							param[ii*7+5], xerror[ii*7+5], 
							param[ii*7+6], xerror[ii*7+6]
							);
					} else {
						printf("\tLine at %f deemed invalid\n", param[ii*7]);
					}
				}
				fflush(fpout);
			} else {
				if (!silent) printf("\tFit deemed invalid, nothing printed to file\n");
			}

			/* Free some memory for next k */
			free(param);
			free(bounds);
			free(xerror);
		}
	}
	fclose(fpout);

	return EXIT_SUCCESS;
}

double model (int numridges, double nu, double k, double theta, double* p)
{
	int ii;
	double ret, den;

	ret = 0.0;
	for (ii=0; ii<numridges; ii++)
	{
		den = nu + k*(p[ii*7+3]*cos(theta)+p[ii*7+4]*sin(theta))/6.28318531 - p[ii*7];
		den = den*den + p[ii*7+2]*p[ii*7+2]/4.0;
		ret += p[ii*7+1]*p[ii*7+2]/(2.0*den);
	}
	ret += p[numridges*7]/(1.0 + pow(nu*p[numridges*7+1], p[numridges*7+2]));
	ret += p[numridges*7+3]*p[numridges*7+5]/(pow(nu-p[numridges*7+4],2.) + pow(0.5*p[numridges*7+5],2.));
/*	ret += p[numridges*5+3]*p[numridges*5+5]/((nu-p[numridges*5+4])*(nu-p[numridges*5+4])+p[numridges*5+5]*p[numridges*5+5]/4.);
*/	return ret;
}

int funk (int m, int n, double* p, double *deviates, double **derivs, void *private)
{
	struct kslice *sub;
	int istw, iendw, num, iw, itht, ik;
	double akt, w1, twopi, tht, den, x2, err, back;

	sub = (struct kslice*) private;
	istw = sub->start;
	iendw = sub->end;
	akt = sub->k;
	twopi = 6.28318530717958647692528676655901; /* lolwut */

	num = 0;
	for (iw=istw; iw<=iendw; iw++)
	{
		w1 = iw*sub->delta_nu;
		back = p[n-3]*p[n-1]/((w1-p[n-2])*(w1-p[n-2]) + 0.25*p[n-1]*p[n-1]);
		back += p[n-6]/(1.+pow(p[n-5]*w1,p[n-4]));
		for (itht=0; itht<sub->ntheta; itht++)
		{
			tht = twopi*itht/sub->ntheta;

			deviates[num] = 0.0;
			if (sub->data[iw-istw][itht] > 0.0 && !isnan(sub->data[iw-istw][itht]))
			{
				for (ik=0; ik<n/7; ik++)
				{
					den = w1 + akt*(p[ik*7+3]*cos(tht)+p[ik*7+4]*sin(tht))/twopi - p[ik*7];
					den = den*den + p[ik*7+2]*p[ik*7+2]/4.0;
					x2 = (p[ik*7+1]+p[ik*7+5]*cos(tht-p[ik*7+6]))*p[ik*7+2]/(2.0*den);
					/*x2 += p[ik*7+5]*cos(tht-p[ik*7+6]);*/
					deviates[num] += x2;
				}
/*				deviates[num] = ((deviates[num] + back) - sub->data[iw-istw][itht]);
				deviates[num] /= sub->noise[iw-istw][itht];
*/
				deviates[num] = sub->data[iw-istw][itht]/(deviates[num]+back);
				deviates[num] -= log(deviates[num]);

			}
			num++;
		}
	}
	return 0;
}

