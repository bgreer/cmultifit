#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include "header.h"

void compute_noise (float*** pol, float*** noise, int nnu, int nk, int ntheta, int radius)
{
	int ii, ij, ik, il;
	float sum, ssum, weight;
	float *gauss;

	/* precompute gaussian kernel */
	gauss = (float*) malloc((6*radius+1)*sizeof(float));
	for (ii=0; ii<6*radius+1; ii++)
		gauss[ii] = exp(-(ii-radius*3)*(ii-radius*3)/(2.*radius*radius));

	for (ij=0; ij<nk; ij++)
	{
		for (ii=0; ii<nnu; ii++)
		{
			sum = weight = 0.0;
			for (il=-radius*3; il<=radius*3; il++)
			{
				if (il+ii>=0 && ii+il<nnu)
				{
					ssum = 0.0;
					for (ik=0; ik<ntheta; ik++)
					{
						ssum += pol[ii+il][ij][ik];
					}
					sum += gauss[il+radius*3]*ssum/ntheta;
					weight += gauss[il+radius*3];
				}
			}
			for (ik=0; ik<ntheta; ik++)
				noise[ii][ij][ik] = 0.0*sum/weight+0.0001;
		}
	}
}

void compute_noise_wavelet (float*** pol, float*** noise, int nnu, int nk, int ntheta, int radius)
{
	int ii, ij, ik, il;
	double* data;
	unsigned int n;
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	n = powerof2(nnu);
	data = (double*) calloc(n,sizeof(double));

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 6);
	work = gsl_wavelet_workspace_alloc (n);
	
	for (ij=0; ij<nk; ij++)
	{
		for (ik=0; ik<ntheta; ik++)
		{
			for (ii=0; ii<nnu; ii++)
				data[ii] = pol[ii][ij][ik];
			for (ii=nnu; ii<n; ii++)
				data[ii] = pol[nnu-1][ij][ik];

			gsl_wavelet_transform_forward (w, data, 1, n, work);
			
			for (ii=0; ii<n/4; ii++)
				data[ii] = 0.0;
			
			gsl_wavelet_transform_inverse (w, data, 1, n, work);
			/* compute stdev */
			for (ii=5; ii<nnu-5; ii++)
			{
				noise[ii][ij][ik] = 0.0;
				for (il=-5; il<=5; il++)
					noise[ii][ij][ik] += (data[ii+il])*(data[ii+il]);
				noise[ii][ij][ik] = sqrt(noise[ii][ij][ik]/11.);
			}
			/* fill ends */
			for (ii=0; ii<5; ii++)
				noise[ii][ij][ik] = noise[5][ij][ik];
			for (ii=nnu-5; ii<nnu; ii++)
				noise[ii][ij][ik] = noise[nnu-6][ij][ik];
		}
	}

	free(data);
}

unsigned int powerof2 (unsigned int n)
{
	int m = 0;
	n--;
	while (n != 1)
	{
		n = n >> 1;
		m++;
	}
	return n << m+1;
}
