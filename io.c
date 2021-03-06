#include <string.h>
#include <math.h>
#include "header.h"

void output_debug (struct params* p, double ***pol, int ntheta, int nk, int nnu, int k, int m, int n, double* x, double delta_nu, double delta_k)
{
	FILE* fp;
	int ii, ik, ij, num;
	double back1, back2, fit, denom, sp, no;
	
	fp = fopen(p->debugfname, "w");
	for (ii=0; ii<nnu; ii++)
	{
		back1 = back2 = fit = sp = no = 0.0;
		for (ik=0; ik<ntheta; ik++)
		{
			back1 += x[n-8]*(1.0+x[n-5]*cos(2.0*(TWOPI*ik/ntheta - x[n-4])))/(1.0+pow(ii*delta_nu/(x[n-7]),x[n-6]));
			back2 += x[n-3]*x[n-1]/((ii*delta_nu-x[n-2])*(ii*delta_nu-x[n-2]) + 0.25*x[n-1]*x[n-1]);
			for (ij=0; ij<(n-NBACK)/NPEAK; ij++)
			{
				denom = ii*delta_nu
					+ k*delta_k*(x[ij*NPEAK+3]*cos(ik*TWOPI/ntheta)
						+x[ij*NPEAK+4]*sin(ik*TWOPI/ntheta))/TWOPI
					- x[ij*NPEAK];
				fit += 0.5*x[ij*NPEAK+2]*x[ij*NPEAK+1] / 
					(denom*denom + 0.25*x[ij*NPEAK+2]*x[ij*NPEAK+2]);
			}
			sp += pol[ii][k][ik];
		}

		fprintf(fp, "%d\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n", ik, ii*delta_nu, 
			sp, 
			0.0, 
			back1+back2+fit, 
			back1, 
			back2, 
			((back1+back2+fit)-sp)/(back1+back2+fit));
	}
	fclose(fp);
}



/* Prints square matrix to output file (Useful for covariance matrix) */
void output_matrix (double* covar, int n, struct params* p)
{
	FILE *fp;
	int ii, ik;
	
	fp = fopen(p->covarfname, "w");
	if (fp==NULL)
	{
		printf("ERROR: Could not open covar output file %s\n", p->covarfname);
		return;
	}
	for (ii=0; ii<n; ii++)
	{
		for (ik=0; ik<n; ik++)
		{
			fprintf(fp, "%d\t%d\t%e\n", ii, ik, covar[ii*n+ik]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}



/* Read FITS spectrum, load header keys into variables, load data cube into spec */
int read_fits_file (double ****spec, struct params* p, 
		int* ntheta, int* nk, int* nnu, double* delta_k, double* delta_nu)
{
	fitsfile *fptr;
	int ii, ij, ik;
	int status = 0;
	long coords[3];
	float *buff, read;

	if (!p->silent) printf("Reading FITS file: %s\n", p->fitsfname);
	
	fits_open_file(&fptr, p->fitsfname, READONLY, &status);
	
	if (status)
	{
		printf("Error opening FITS file: %s\n", p->fitsfname);
		fits_report_error(stdout, status);
		return EXIT_FAILURE;
	}
	
	fits_read_key(fptr, TINT, "NAXIS1", ntheta, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS2", nk, NULL, &status);
	fits_read_key(fptr, TINT, "NAXIS3", nnu, NULL, &status);
	
	if (status)
	{
		printf("Could not find FITS dimensions.\n");
		printf("Make sure the following keys are set: NAXIS1, NAXIS2, NAXIS3\n");
		return EXIT_FAILURE;
	}
	
	if (!p->silent)
		printf("FITS dimensions:\n\ttheta:\t%d\n\tk:\t%d\n\tnu:\t%d\n", *ntheta, *nk, *nnu);
	
	fits_read_key(fptr, TFLOAT, "DELTA_K", &read, NULL, &status);
	*delta_k = read;
	fits_read_key(fptr, TFLOAT, "DELTA_NU", &read, NULL, &status);
	*delta_nu = read;
	
	if (status)
	{
		printf("Could not find FITS keys DELTA_NU and/or DELTA_K.\n");
		return EXIT_FAILURE;
	}
	if (!p->silent) printf("\tDELTA_NU:\t%f\n\tDELTA_K:\t%f\n", *delta_nu, *delta_k);

/* Allocate memory for entire FITS data cube */
	if (!p->silent) 
		printf("\nAllocating %ld bytes for data cube.\n", sizeof(double)*(*ntheta)*(*nk)*(*nnu));
	if (!p->silent)
		printf("\t With %ld bytes of overhead.\n", sizeof(double*)*(*nk)*(*nnu)+sizeof(double**)*(*nnu));

	*spec = (double***) malloc((*nnu)*sizeof(double**));
	for (ii=0; ii<(*nnu); ii++)
	{
		(*spec)[ii] = (double**) malloc((*nk)*sizeof(double*));
		for (ij=0; ij<(*nk); ij++)
		{
			(*spec)[ii][ij] = (double*) malloc((*ntheta)*sizeof(double));
			if ((*spec)[ii][ij]==NULL)
			{
				printf("Error allocating memory.\n");
				return EXIT_FAILURE;
			}
		}
	}
	buff = (float*) malloc((*ntheta)*sizeof(float));
/* Read FITS file into array */
	coords[0] = 1L;
	if (!p->silent) printf("Reading data cube into memory.\n");
	for (coords[2]=1; coords[2]<=(*nnu); coords[2]++)
	{
		for (coords[1]=1; coords[1]<=(*nk); coords[1]++)
		{
			fits_read_pix(fptr, TFLOAT, coords, (*ntheta), NULL, buff, NULL, &status);
			for (ik=0; ik<(*ntheta); ik++)
			{
				if (buff[ik] < 0.0 || isnan(buff[ik])) buff[ik] = 0.0;
				(*spec)[coords[2]-1][coords[1]-1][ik] = buff[ik];
			}
		}
	}
	fits_close_file(fptr, &status);
	free(buff);
}

/* Read model file */
int read_model_file (struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width)
{
	FILE *fpmodel;
	int ii, ij, ik, il;
	float readfreq, readamp, readwidth, readk;

	if (!par->silent) printf("Reading model from %s\n", par->modelfname);
	*numridges = (int*) malloc(nk*sizeof(int));
	*freq = (double**) malloc(nk*sizeof(double*));
	*amp = (double**) malloc(nk*sizeof(double*));
	*width = (double**) malloc(nk*sizeof(double*));
	for (ii=0; ii<nk; ii++)
		(*numridges)[ii] = 0;
	/* Run through once to count ridges */
	fpmodel = fopen(par->modelfname, "r");
	if (fpmodel==NULL)
	{
		printf("ERROR: could not open model file: %s\n", par->modelfname);
		return EXIT_FAILURE;
	}
	while (fscanf(fpmodel, "%d\t%e\f%e\t%e\t%e\n", &ii, &readk, &readfreq, &readamp, &readwidth) != EOF)
	{
		if (ii>=0 && ii<nk)
		{
			(*numridges)[ii]++;
		} else {
			printf("ERROR: k-value out of range in model: %d\n", ii);
			return EXIT_FAILURE;
		}
	}
	rewind(fpmodel);
	/* allocate enough for all ridges, then load them */
	for (ii=0; ii<nk; ii++)
	{
		if ((*numridges)[ii]>0)
		{
			(*freq)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*amp)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
			(*width)[ii] = (double*) malloc((*numridges)[ii]*sizeof(double));
		}
		for (ij=0; ij<(*numridges)[ii]; ij++)
		{
			il = fscanf(fpmodel, "%d\t%e\t%e\t%e\t%e\n", &ik, &readk, &readfreq, &readamp, &readwidth);
			if (readfreq<0.0 || readamp<0.0 || readwidth<0.0)
			{
				printf("ERROR: Invalid parameters in model, k=%d\n", ii);
				return EXIT_FAILURE;
			}
			(*freq)[ii][ij] = readfreq;
			(*amp)[ii][ij] = readamp;
			(*width)[ii][ij] = readwidth*0.7;
		}
	}
	fclose(fpmodel);
}


/* Reads parameter file fname and loads data into parameter struct *p */
void read_param_file (char* fname, struct params* p)
{
	FILE *fp;
	char buffer[200];
	int index, ii;

	fp = fopen(fname, "r");
	if (fp==NULL)
	{
		printf("Could not open parameter file: %s\n", fname);
		return;
	}

	/* Allocate space for filenames */
	p->fitsfname = malloc(200*sizeof(char));
	p->modelfname = malloc(200*sizeof(char));
	p->outfname = malloc(200*sizeof(char));
	p->debugfname = p->covarfname = 0;

	/* Scan through for non-comment and non-blank lines */
	index = 0;
	while (fgets(buffer, sizeof(buffer), fp) != NULL)
	{
		if (buffer[0] != '#' && buffer[0] != '\n')
		{
			trim(buffer);
			switch (index)
			{
				case PARAM_SPECTRUM:
					strcpy(p->fitsfname, buffer);
					break;
				case PARAM_MODEL:
					strcpy(p->modelfname, buffer);
					break;
				case PARAM_OUTPUT:
					strcpy(p->outfname, buffer);
					break;
				case PARAM_SILENT:
					p->silent = atoi(buffer);
					break;
				case PARAM_KRANGE:
					/* split string by putting null-terminator where the space is */
					ii=0;
					while (ii<strlen(buffer) && buffer[ii]!=' ')
						ii++;
					buffer[ii] = 0;
					p->kstart = atoi(buffer);
					p->kend = atoi(buffer+ii+1); /* pointer arithmetic */
					break;
				case PARAM_FTOL:
					p->ftol = atof(buffer);
					break;
				case PARAM_XTOL:
					p->xtol = atof(buffer);
					break;
				case PARAM_GTOL:
					p->gtol = atof(buffer);
					break;
				case PARAM_NITER:
					p->niter = atoi(buffer);
					break;
				case PARAM_DEBUG:
					if (strcmp(buffer, "0")!=0)
					{
						p->debugfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->debugfname, buffer);
					}
					break;
				case PARAM_COVAR:
					if (strcmp(buffer, "0")!=0)
					{
						p->covarfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->covarfname, buffer);
					}
					break;
				case PARAM_BACK:
					if (strcmp(buffer, "0")!=0)
					{
						p->backfname = malloc(strlen(buffer)*sizeof(char));
						strcpy(p->backfname, buffer);
					}
					break;
				case PARAM_FIT:
					p->dofits = atoi(buffer);
					break;
				case PARAM_CUTOFF:
					p->ac_cutoff = (double) atof(buffer);
					break;
				case PARAM_FITABOVE:
					p->fit_above = atoi(buffer);
					break;
				case PARAM_RIDGE:
					p->detect_ridge = atoi(buffer);
					break;
			}
			index++;
		}
	}
	fclose(fp);

	/* If not silent, print run params (for sanity checks) */
	if (!p->silent)
	{
		printf("\nRun parameters:\n");
		printf("\tAcoustic Cutoff Frequency: %6.1f uHz\n", p->ac_cutoff);
		switch (p->fit_above)
		{
			case 0:
				printf("\tNo peaks fit above this.\n"); break;
			case 1:
				printf("\tSimple peaks fit above this.\n"); break;
			case 2:
				printf("\tFull peaks fit above this.\n"); break;
		}
		if (p->debugfname)
			printf("\tDebug filename base: %s\n", p->debugfname);
		if (p->covarfname)
			printf("\tCovar filename base: %s\n", p->covarfname);
		printf("\n");
	}
}



/* fgets leaves a trailing newline, this removes it */
void trim (char* str)
{
	int len;
	len = strlen(str)-1;
	if (str[len] == '\n')
		str[len] = 0;
}
