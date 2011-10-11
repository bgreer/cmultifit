#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "header.h"


/* Computes gaussian smoothed spectrum in nu (azimuthally averaged) 
	sets noise to be 0.1*(local avg)
	radius = sigma of gaussian smoothing kernel */
void compute_noise_smooth (double*** pol, double*** noise, 
						int nnu, int nk, int ntheta, int radius)
{
	int ii, ij, ik, il;
	double sum, ssum, weight;
	double *gauss;

	/* precompute gaussian kernel */
	gauss = (double*) malloc((6*radius+1)*sizeof(double));
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
				noise[ii][ij][ik] = 0.1*sum/weight;
		}
	}
	free(gauss);
}

/* Computes noise to be 0.1*(average power at constant k) */
void compute_noise_const (double*** pol, double*** noise, int nnu, int nk, int ntheta)
{
	int ii, ij, ik;
	double sum;
	for (ij=0; ij<nk; ij++)
	{
		sum = 0.0;
		for (ii=0; ii<nnu; ii++)
		{
			for (ik=0; ik<ntheta; ik++)
				sum += pol[ii][ij][ik];
		}
		sum /= 10.*(ntheta*nnu);
		for (ii=0; ii<nnu; ii++)
			for (ik=0; ik<ntheta; ik++)
				noise[ii][ij][ik] = sum;
	}
}


/* Wavelet transform requires data to be 2^m length
	this determines next highest 2^m >= n
	bit-shifting makes code look more legit than it is */
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
