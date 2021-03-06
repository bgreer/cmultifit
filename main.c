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
	FILE *fpout, *fpback;
	struct params par;
	int status = 0;
	int ii, ij, ik;
	int ntheta, nk, nnu;
	double delta_nu, delta_k;
	double ***pol;
	double **freq, **amp, **width;
	int *numridges, **fit_type, *ridgenum;

	int mpreturn, index;
	struct kslice subsection;
	double* param;
	mp_result *mpres;
	mp_config *mpconf;
	mp_par *bounds;
	double *xerror;
	double *covar;
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
	if (read_fits_file(&pol, &par, &ntheta, &nk, &nnu, &delta_k, &delta_nu)==EXIT_FAILURE)
		return EXIT_FAILURE;
	thtarr = malloc(ntheta*sizeof(double));
	thtpow = malloc(ntheta*sizeof(double));

	/* Check kstart and kend */
	if (par.kstart > par.kend) par.kend = par.kstart;
	if (par.kstart < 0) par.kstart = 0;
	if (par.kstart >= nk) par.kstart = nk-1;
	if (par.kstart < 0) par.kend = 0;
	if (par.kend >= nk) par.kend = nk-1;

	/* Read in model file verbatim */
	read_model_file(&par, nk, &numridges, &freq, &amp, &width);

	/* Decide which peaks to fit in what way */
	fit_type = malloc(nk*sizeof(int*));
	for (ii=0; ii<nk; ii++)
	{
		fit_type[ii] = 0;
		if (numridges[ii] > 0)
		{
			fit_type[ii] = malloc(numridges[ii]*sizeof(int));
			for (ij=0; ij<numridges[ii]; ij++)
			{
				if (par.fit_above==2)
					fit_type[ii][ij] = 1; /* Full line profile */
				else if (par.fit_above==1 && freq[ii][ij] >= par.ac_cutoff)
					fit_type[ii][ij] = 2; /* Simple line profile */
				else if (par.fit_above==0 && freq[ii][ij] >= par.ac_cutoff)
					fit_type[ii][ij] = 0; /* Do not fit */
				else
					fit_type[ii][ij] = 1;
			}
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


/* BEGIN LOOP OVER K-BINS */

	for (ij=par.kstart; ij<=par.kend; ij++)
	{
		if (numridges[ij]>0)
		{
			if (!par.silent) printf("Starting fit of %d ridges at k=%d\n", numridges[ij], ij);
			param = (double*) malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			bounds = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(mp_par));
			xerror = malloc((numridges[ij]*NPEAK+NBACK)*sizeof(double));
			covar = malloc((numridges[ij]*NPEAK+NBACK)*(numridges[ij]*NPEAK+NBACK)*sizeof(double*));
			ridgenum = malloc(numridges[ij]*sizeof(int));

			/* Do rough fit of single peaks */
			if (!par.silent) printf("\tDoing single ridge estimates.\n");
			for (ii=0; ii<numridges[ij]; ii++)
			{
				if (par.dofits)
					fit_peak(&par, &(freq[ij][ii]), &(amp[ij][ii]), &(width[ij][ii]), pol, 
						delta_nu, delta_k, nnu, ntheta, ij);
			}
			
			/* Do rough fit of background */
			if (!par.silent) printf("\tFitting background at low frequency.\n");
			fit_back(&par, &(param[numridges[ij]*NPEAK]), &(param[numridges[ij]*NPEAK+1]), 
				&(param[numridges[ij]*NPEAK+2]), pol, delta_nu, nnu, ntheta, ij);

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
					bounds[ii*NPEAK+ik].deriv_reltol = 0.0;
					bounds[ii*NPEAK+ik].deriv_abstol = 0.0;
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
				if (param[ii*NPEAK] < par.ac_cutoff || par.fit_above != 1)
				{
					param[ii*NPEAK+3] = 1.0;
					bounds[ii*NPEAK+3].limited[0] = bounds[ii*NPEAK+3].limited[1] = 1;
					bounds[ii*NPEAK+3].limits[0] = -1000.0;
					bounds[ii*NPEAK+3].limits[1] = 1000.0;
					param[ii*NPEAK+4] = 1.0;
					bounds[ii*NPEAK+4].limited[0] = bounds[ii*NPEAK+4].limited[1] = 1;
					bounds[ii*NPEAK+4].limits[0] = -1000.0;
					bounds[ii*NPEAK+4].limits[1] = 1000.0;
				} else {
					param[ii*NPEAK+3] = 0.;
					param[ii*NPEAK+4] = 0.;
					bounds[ii*NPEAK+3].fixed = 1;
					bounds[ii*NPEAK+4].fixed = 1;
					bounds[ii*NPEAK+3].limited[0] = bounds[ii*NPEAK+3].limited[1] = 0;
					bounds[ii*NPEAK+4].limited[0] = bounds[ii*NPEAK+4].limited[1] = 0;
				}
				/* set theta variation params */
				param[ii*NPEAK+5] = 0.01;
				bounds[ii*NPEAK+5].limited[0] = bounds[ii*NPEAK+5].limited[1] = 1;
				bounds[ii*NPEAK+5].limits[0] = -1.0;
				bounds[ii*NPEAK+5].limits[1] = 1.0;
				param[ii*NPEAK+6] = PI/2.0;
				bounds[ii*NPEAK+6].limited[0] = bounds[ii*NPEAK+6].limited[1] = 0;

				for (ik=0; ik<7; ik++)
					if (bounds[ii*NPEAK+ik].limited[0] && param[ii*NPEAK+ik] < bounds[ii*NPEAK+ik].limits[0])
						printf("error on param %d of peak %d\t%e %e\n", ik, ii, param[ii*NPEAK+ik], bounds[ii*NPEAK+ik].limits[0]);
			}

			for (ii=0; ii<NBACK; ii++)
			{
				bounds[numridges[ij]*NPEAK+ii].fixed = 0;
				bounds[numridges[ij]*NPEAK+ii].side = 0;
				bounds[numridges[ij]*NPEAK+ii].deriv_debug = 0;
				bounds[numridges[ij]*NPEAK+ii].relstep = 0.0;
				bounds[numridges[ij]*NPEAK+ii].step = 0.0;
				bounds[numridges[ij]*NPEAK+ii].deriv_reltol = 0.0;
				bounds[numridges[ij]*NPEAK+ii].deriv_abstol = 0.0;
			}

			/* Values and bounds for background power law */
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

			/* Values and bounds for background lorentzian */
			param[numridges[ij]*NPEAK+5] = 0.03*param[1];
			bounds[numridges[ij]*NPEAK+5].limited[0] = 1;
			bounds[numridges[ij]*NPEAK+5].limited[1] = 0;
			bounds[numridges[ij]*NPEAK+5].limits[0] = 0.0;
			param[numridges[ij]*NPEAK+6] = 1500.;
			bounds[numridges[ij]*NPEAK+6].limited[0] = bounds[numridges[ij]*NPEAK+6].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+6].limits[0] = 200.;
			bounds[numridges[ij]*NPEAK+6].limits[1] = 4000;/*nnu*delta_nu;*/
			param[numridges[ij]*NPEAK+7] = 1000.;
			bounds[numridges[ij]*NPEAK+7].limited[0] = bounds[numridges[ij]*NPEAK+7].limited[1] = 1;
			bounds[numridges[ij]*NPEAK+7].limits[0] = 500.0;
			bounds[numridges[ij]*NPEAK+7].limits[1] = 5000.;

			/* for brad's fake spectrum, fix background to 0 */
			/*
			param[numridges[ij]*NPEAK] = 0.0;
			bounds[numridges[ij]*NPEAK].fixed = 1;
			param[numridges[ij]*NPEAK+5] = 0.0;
			bounds[numridges[ij]*NPEAK+5].fixed = 1;
			*/

			subsection.par = &par;
			subsection.n = numridges[ij];
			
			/* Load klice struct for passing to function */
			subsection.start = 100/delta_nu;
			/* testing? */
			/*
			subsection.start = (param[0]-param[2])/delta_nu;
			if (subsection.start < 0) subsection.start = 0;
			*/
			subsection.end = (param[(numridges[ij]-1)*NPEAK]+1.*param[(numridges[ij]-1)*NPEAK+2])/delta_nu;
			if (subsection.end > nnu-1) subsection.end = nnu-1;

			subsection.delta_nu = delta_nu;
			subsection.ntheta = ntheta;
			subsection.k = (ij+1)*delta_k;
			subsection.data = malloc((subsection.end-subsection.start+1)*sizeof(double*));
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				subsection.data[ii-subsection.start] = malloc(ntheta*sizeof(double));
				for (ik=0; ik<ntheta; ik++)
				{
					subsection.data[ii-subsection.start][ik] = pol[ii][ij][ik];
				}
			}

			/* Set optimization parameters */
			mpconf->ftol = par.ftol;
			mpconf->xtol = par.xtol;
			mpconf->gtol = par.gtol;
			mpconf->covtol = 1e-14;
			mpconf->maxiter = par.niter;
			mpconf->nofinitecheck = 1;
			mpreturn = 1;

/* Perform optimization */
			if (par.dofits && !par.silent) printf("\tDoing optimization...\n");			
			if (par.dofits) mpreturn = mpfit(&funk, 
					ntheta*(subsection.end-subsection.start+1), 
					numridges[ij]*NPEAK+NBACK, 
					param, 
					bounds, 
					mpconf, 
					&subsection, 
					mpres);


			/* Check return value of mpfit */
			mpreturn_translate(mpreturn);
		
			/* Output optimization details */
			if (!par.silent) 
			{
				printf("\tStarting LH = %e\n", 
					mpres->orignorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tFinal LH = %e\n", 
					mpres->bestnorm/(ntheta*(subsection.end-subsection.start)));
				printf("\tNumber of iterations = %d\n", mpres->niter);
				printf("\tNumber of function evals = %d\n", mpres->nfev);
				printf("\tNumber of pegged parameters = %d\n", mpres->npegged);
			}


			/* Determine which parameters are pegged */
			if (mpres->npegged)
			{
				for (ii=0; ii<numridges[ij]; ii++)
				{
					for (ik=0; ik<NPEAK; ik++)
					{
						if (bounds[ii*NPEAK+ik].limited[0] && param[ii*NPEAK+ik] == bounds[ii*NPEAK+ik].limits[0])
							printf("Param %d of peak %d pegged low\n", ik, ii);
						if (bounds[ii*NPEAK+ik].limited[1] && param[ii*NPEAK+ik] == bounds[ii*NPEAK+ik].limits[1])
							printf("Param %d of peak %d pegged high\n", ik, ii);
					}
				}
				for (ik=0; ik<NBACK; ik++)
				{
					if (bounds[numridges[ij]*NPEAK+ik].limited[0] && 
							param[numridges[ij]*NPEAK+ik] == bounds[numridges[ij]*NPEAK+ik].limits[0])
						printf("Param %d of back pegged low\n", ik);
					if (bounds[numridges[ij]*NPEAK+ik].limited[1] && 
							param[numridges[ij]*NPEAK+ik] == bounds[numridges[ij]*NPEAK+ik].limits[1])
						printf("Param %d of back pegged high\n", ik);
				}
			}


			/* Compute error bars */
			errbars(numridges[ij]*NPEAK+NBACK, param, &subsection, &covar);
/*			for (ii=0; ii<numridges[ij]*NPEAK+NBACK; ii++)
			{
				for (ik=0; ik<numridges[ij]*NPEAK+NBACK; ik++)
				{
					printf("%d\t%d\t%e\n", ii, ik, covar[ii*(numridges[ij]*NPEAK+NBACK)+ik] / 
							(param[ii]*param[ik]));
				}
				printf("\n");
			}
*/			
			/* Load covar into xerror because I'm lazy */
			for (ii=0; ii<numridges[ij]*NPEAK+NBACK; ii++)
				xerror[ii] = sqrt(fabs(covar[ii*(numridges[ij]*NPEAK+NBACK)+ii]));

			/* Post-process anisotropy parameters (negative amplitude -> pi/2 shift) */
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

			/* Auto-detect ridge */
			if (par.detect_ridge)
			{
				for (ii=0; ii<numridges[ij]; ii++)
					ridgenum[ii] = detect_ridge((ij+1)*delta_k, param[ii*NPEAK]);
			} else {
				for (ii=0; ii<numridges[ij]; ii++)
					ridgenum[ii] = 0;
			}

			/* Print fit debug */
			if (par.debugfname) output_debug(&par, pol, ntheta, nk, nnu, ij, ntheta*(subsection.end-subsection.start+1), 
					numridges[ij]*NPEAK+NBACK, 
					param, delta_nu, delta_k);


			/* Print covariance matrix */
/*			if (par.covarfname) output_matrix(covar, numridges[ij]*NPEAK+NBACK, &par);*/


			/* Output fit to file if valid */
			if (mpreturn!=MP_MAXITER && mpreturn > 0) /* Removed 'invalid if mpreturn=3' */
			{
				for (ii=0; ii<numridges[ij]; ii++)
				{
					if (param[ii*NPEAK] <= par.ac_cutoff || 1)
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
					}/* else if (xerror[ii*NPEAK]<=0.0 || xerror[ii*NPEAK+1]<=0.0 || 
							xerror[ii*NPEAK+2]<=0.0 || xerror[ii*NPEAK+3]<=0.0 || 
							xerror[ii*NPEAK+4]<=0.0) {
						printf("\tLine at %f deemed invalid: ", param[ii*NPEAK]);
						printf("Zero error\n");
					}*/ else {
						fprintf(fpout, "%d\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%f\t%e\t%e\t%e\t%d\n", 
							ij, 
							(ij+1)*delta_k,
							param[ii*NPEAK]/((ij+1)*delta_k), 
							param[ii*NPEAK], xerror[ii*NPEAK+3], 
							param[ii*NPEAK+1], xerror[ii*NPEAK+1], 
							param[ii*NPEAK+2], xerror[ii*NPEAK+2], 
							param[ii*NPEAK+3], xerror[ii*NPEAK+3], 
							param[ii*NPEAK+4], xerror[ii*NPEAK+4],
							param[ii*NPEAK+5], xerror[ii*NPEAK+5], 
							param[ii*NPEAK+6], xerror[ii*NPEAK+6], 
							param[ii*NPEAK+7], xerror[ii*NPEAK+7],
							ridgenum[ii]);
					}
					}
				}
				fflush(fpout);
			

				/* Output background params to separate file */
				if (par.backfname)
				{
					fprintf(fpback, "%d\t%f", ij, ij*delta_k);
					for (ik=0; ik<NBACK; ik++)
						fprintf(fpback, "\t%e\t%e", param[numridges[ij]*NPEAK+ik], 
							xerror[numridges[ij]*NPEAK+ik]);
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
			for (ii=subsection.start; ii<=subsection.end; ii++)
			{
				free(subsection.data[ii-subsection.start]);
			}
			free(subsection.data);
			free(covar);
			free(ridgenum);
		}
	}

/* END LOOP OVER K-BINS */


	fclose(fpout);
	if (par.backfname) fclose(fpback);
	
	/* Free ALL THE THINGS */
	for (ii=0; ii<nk; ii++)
	{
		if (numridges[ii]>0)
		{
			free(freq[ii]);
			free(amp[ii]);
			free(width[ii]);
		}
	}
	free(freq);
	free(amp);
	free(width);
	free(numridges);
	free(mpres);
	free(mpconf);
	free(thtarr);
	free(thtpow);
	free(par.fitsfname);
	free(par.modelfname);
	free(par.outfname);
	if (par.debugfname) free(par.debugfname);
	if (par.covarfname) free(par.covarfname);
	if (par.backfname) free(par.backfname);
	for (ii=0; ii<nnu; ii++)
	{
		for (ij=0; ij<nk; ij++)
		{
			free(pol[ii][ij]);
		}
		free(pol[ii]);
	}
	free(pol);
	return EXIT_SUCCESS;
}

void mpreturn_translate (int mpreturn)
{
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
}
