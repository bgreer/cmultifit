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
				case 0: /* FITS spectrum */
					strcpy(p->fitsfname, buffer);
					break;
				case 1: /* model file */
					strcpy(p->modelfname, buffer);
					break;
				case 2: /* output file */
					strcpy(p->outfname, buffer);
					break;
				case 3: /* silent mode */
					p->silent = atoi(buffer);
					break;
			}
			index++;
		}
	}
	fclose(fp);
}

/* fgets leaves a trailing newline, this removes it */
void trim (char* str)
{
	int len;
	len = strlen(str)-1;
	if (str[len] == '\n')
		str[len] = 0;
}
