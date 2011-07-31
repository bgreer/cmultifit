#include <string.h>
#include "header.h"

/* Reads parameter file fname and loads data into parameter struct *p */
void read_param_file (char* fname, struct params* p)
{
	FILE *fp;
	char buffer[200];
	int index;

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
