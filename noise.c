#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
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

/* Computes noise to be stdev of reconstructed spectrum at constant (k, theta)
	using only smallest scale wavelet components*/
void compute_noise_wavelet (double*** pol, double*** noise, int nnu, int nk, int ntheta)
{
	int ii, ij, ik, il, offset;
	double* data, frontavg, backavg, count;
	unsigned int n;
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	n = powerof2(nnu);
	offset = (n-nnu)/2;
	data = (double*) calloc(n,sizeof(double));
	
	/* prepare for wavelet transforms */
	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (n);
	for (ij=0; ij<nk; ij++)
	{
		for (ik=0; ik<ntheta; ik++)
		{
			/* find frontavg, backavg */
			frontavg = backavg = 0.0;
			count = 0.0;
			for (ii=0; ii<10; ii++)
			{
				if (pol[ii][ij][ik] > 0.0)
				{
					count += 1.0;
					frontavg += pol[ii][ij][ik];
				}
			}
			frontavg /= count;
			count = 0.0;
			for (ii=nnu-10; ii<nnu; ii++)
			{
				if (pol[ii][ij][ik] > 0.0)
				{
					count += 1.0;
					backavg += pol[ii][ij][ik];
				}
			}
			backavg /= count;
			/* fill data array with spectrum data, pad ends */
			for (ii=0; ii<nnu; ii++)
				data[ii+offset] = pol[ii][ij][ik];
			for (ii=0; ii<offset; ii++)
				data[ii] = frontavg;
			for (ii=nnu+offset; ii<n; ii++)
				data[ii] = backavg;

			for (ii=0; ii<nnu; ii++)
			{
				if (data[ii] == 0.0) data[ii] = data[ii+1];
			}

			/* do wavelet transform */
			gsl_wavelet_transform_forward (w, data, 1, n, work);
			
			/* filter in wavelet space(?) */
			for (ii=0; ii<n/4; ii++)
				data[ii] = 0.0;
			
			/* transform back */
			gsl_wavelet_transform_inverse (w, data, 1, n, work);
/*if (ik==0 && ij==20)
	for (ii=0; ii<n; ii++)
		printf("%d\t%e\n", ii, data[ii]);
*/			/* compute stdev */
			for (ii=4; ii<nnu-4; ii++)
			{
				noise[ii][ij][ik] = 0.0;
				for (il=-4; il<=4; il++)
					noise[ii][ij][ik] += (data[ii+offset+il])*(data[ii+offset+il]);
				noise[ii][ij][ik] = sqrt(noise[ii][ij][ik]/9.);
			}
			/* fill ends */
			for (ii=0; ii<4; ii++)
				noise[ii][ij][ik] = noise[4][ij][ik];
			for (ii=nnu-4; ii<nnu; ii++)
				noise[ii][ij][ik] = noise[nnu-5][ij][ik];
		}
	}
	free(data);
	gsl_wavelet_free(w);
	gsl_wavelet_workspace_free(work);
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
