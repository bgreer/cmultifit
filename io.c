#include <string.h>
#include "header.h"

/* Read FITS spectrum, load header keys into variables, load data cube into spec */
int read_fits_file (float**** spec, float**** noise, struct params* p, 
		int* ntheta, int* nk, int* nnu, float* delta_k, float* delta_nu)
{
	fitsfile *fptr;
	int ii, ij, ik;
	int status = 0;
	long coords[3];
	float *buff;

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
	
	fits_read_key(fptr, TFLOAT, "DELTA_K", delta_k, NULL, &status);
	fits_read_key(fptr, TFLOAT, "DELTA_NU", delta_nu, NULL, &status);
	
	if (status)
	{
		printf("Could not find FITS keys DELTA_NU and/or DELTA_K.\n");
		return EXIT_FAILURE;
	}
	if (!p->silent) printf("\tDELTA_NU:\t%f\n\tDELTA_K:\t%f\n", *delta_nu, *delta_k);

/* Allocate memory for entire FITS data cube */
	if (!p->silent) 
		printf("\nAllocating %ld bytes for data cube.\n", sizeof(float)*(*ntheta)*(*nk)*(*nnu));
	if (!p->silent)
		printf("\t With %ld bytes of overhead.\n", sizeof(float*)*(*nk)*(*nnu) + sizeof(float**)*(*nnu));

	*spec = (float***) malloc((*nnu)*sizeof(float**));
	*noise = (float***) malloc((*nnu)*sizeof(float**));
	for (ii=0; ii<(*nnu); ii++)
	{
		(*spec)[ii] = (float**) malloc((*nk)*sizeof(float*));
		(*noise)[ii] = (float**) malloc((*nk)*sizeof(float*));
		for (ij=0; ij<(*nk); ij++)
		{
			(*spec)[ii][ij] = (float*) malloc((*ntheta)*sizeof(float));
			(*noise)[ii][ij] = (float*) malloc((*ntheta)*sizeof(float));
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
			}
			memcpy((*spec)[coords[2]-1][coords[1]-1], buff, (*ntheta)*sizeof(float));
		}
	}
	fits_close_file(fptr, &status);
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

	/* Scan through for non-comment and non-blank lines */
	index = 0;
	while (fgets(buffer, sizeof(buffer), fp) != NULL)
	{
		if (buffer[0] != '#' && buffer[0] != '\n')
		{
			trim(buffer);
			switch (index)
			{
				case PARAM_SPECTRUM: /* FITS spectrum */
					strcpy(p->fitsfname, buffer);
					break;
				case PARAM_MODEL: /* model file */
					strcpy(p->modelfname, buffer);
					break;
				case PARAM_OUTPUT: /* output file */
					strcpy(p->outfname, buffer);
					break;
				case PARAM_SILENT: /* silent mode */
					p->silent = atoi(buffer);
					break;
				case PARAM_WEIGHT: /* weighting */
					p->chiweight = atoi(buffer);
					break;
				case PARAM_NOISE: /* noise mode */
					p->noisemode = atoi(buffer);
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
			}
			index++;
		}
	}
	fclose(fp);

	/* If not silent, print run params (for sanity checks) */
	if (!p->silent)
	{
		printf("\nRun parameters:\n");
		printf("\tChi-squared weighting: ");
		switch (p->chiweight)
		{
			case WEIGHT_NOISE:
				printf("Noise Weighted\n");
				break;
			case WEIGHT_MAXL:
				printf("Fit Weighted (Maximum Likelihood)\n");
				break;
		}
		if (p->chiweight == WEIGHT_NOISE)
		{
			printf("\tNoise estimator: ");
			switch (p->noisemode)
			{
				case NOISE_CONST:
					printf("Constant Noise\n");
					break;
				case NOISE_SMOOTH:
					printf("Smooth Noise\n");
					break;
				case NOISE_WAVELET:
					printf("Wavelet Noise\n");
					break;
			}
		}
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
