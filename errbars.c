#include "header.h"
#include "lapacke.h"
#include <math.h>

/* This is going to suck */


/* 
	INPUT:
		all data for a slice
		final parameters (none pegged)
			(background can be pegged, just ignore it?)

	OUTPUT:
		covariance matrix

	PROCESS:
		I'll tell you when it's done
*/


/* Main function to be called */
void errbars (int numparams, double *p, struct kslice *ks, double **covar)
{
	int ii, ij, ik;
	lapack_int n, lda, info, worksize;
	char c;
	double *fish, *fish2, *work, *rscale, *cscale;
	double row, col, amax;
	int *piv;

	printf("\tBEGIN ERROR BAR COMPUTATION\n\tUSER ENCOURAGED TO PRAY\n");

	c = 'U';
	fish = malloc(numparams*numparams*sizeof(double));
	piv = malloc(numparams*numparams*sizeof(int));
	fish2 = malloc(numparams*numparams*sizeof(double));
	rscale = malloc(numparams*sizeof(double));
	cscale = malloc(numparams*sizeof(double));

	/* compute fisher information matrix */
	fisher(numparams, p, ks, fish);
	fisher(numparams, p, ks, fish2);

	n = numparams;
	lda = numparams;
	dgeequ_(&n, &n, fish, &lda, rscale, cscale, &row, &col, &amax, &info);
/*
	for (ii=0; ii<numparams; ii++)
		for (ij=0; ij<numparams; ij++)
			fish[ii*numparams+ij] *= rscale[ii]*cscale[ij];
*/

	n = numparams;
	lda = numparams;
	dgetrf_(&n, &n, fish, &lda, piv, &info);

	worksize = 32*n;
	work = malloc(worksize*sizeof(double));
	dgetri_(&n, fish, &lda, piv, work, &worksize, &info);
/*
	for (ii=0; ii<numparams; ii++)
		for (ij=0; ij<numparams; ij++)
			fish[ii*numparams+ij] *= rscale[ij]*cscale[ii];
*/
	/* compute inverse of fisher information matrix */
/*
	for (ii=0; ii<numparams; ii++)
	{
		for (ij=0; ij<numparams; ij++)
			printf("%d\t%d\t%e\n", ii, ij, (fish[ii*numparams+ij]));
		printf("\n");
	}
*/
	/* return */
	*covar = fish;

	/* free local memory */
	free(work);
	free(piv);
}


/* Compute the fisher information matrix
	each element is the 2nd derivative of the log-likelihood evaluated at convergence
	element i,j = ( d^2/(dpi*dpj)[ln(f(x,p))] ) | p_f
*/
void fisher (int numparams, double *p, struct kslice *ks, double *A)
{
	int ii, ij, inu, itht;
	double corr1, corr2, mod, datapt;
	double *d1, **d2;


	d1 = calloc(numparams,sizeof(double));
	d2 = malloc(numparams*sizeof(double*));
	for (ii=0; ii<numparams; ii++)
		d2[ii] = calloc(numparams,sizeof(double));

	/* clear matrix in perperation of incremental sums */
	for (ii=0; ii<numparams*numparams; ii++)
		A[ii] = 0.0;

	/* loop over data, summing ML and matrix elements */
	for (inu=ks->start; inu<=ks->end; inu++)
	{
		for (itht=0; itht<ks->ntheta; itht++)
		{
			datapt = ks->data[inu-ks->start][itht];
			/* compute things for a single pixel */
			mod = model(p, ks, inu, itht);
			corr1 = (1.0 - datapt/mod)/mod;
			corr2 = (2.0*datapt/mod - 1.0)/(mod*mod);
			
			/* call to find all first derivatives */
			dm(p, ks, inu, itht, d1);
			d2m(p, ks, inu, itht, d2);
			
			/* matrix element is d2m/didj*corr1 + dm/di*dm/dj*corr2 */
			for (ii=0; ii<numparams; ii++)
			{
				for (ij=0; ij<numparams; ij++)
				{
					/* add contribution to matrix element ii,ij */
					A[ii*numparams+ij] += d2[ii][ij]*corr1 + d1[ii]*d1[ij]*corr2;
				}
			}
		}
	}

	for (ii=0; ii<numparams; ii++)
	{
		for (ij=0; ij<numparams; ij++)
		{
			/* load data */
			if (A[ii*numparams+ij] != A[ij*numparams+ii])
				printf("problem at %d %d\n", ii, ij);
			/*A[ij*numparams+ii] = A[ij*numparams+ii];*/
		}
	}
	free(d1);
}

/* compute all cross derivatives at a single point */
/* there is a clever way of doing this TODO: find it */
void d2m (double *p, struct kslice *ks, int inu, int itht, double **d)
{
	int ii, ij, nr, offset;
	double lor, w1, tht, akt, back, backw, back2, den, denp, co, twot, anis, shift, den1, A, G;
	double po;

	nr = ks->n;
	offset = nr*NPEAK;

	for (ii=0; ii<nr*NPEAK+NBACK; ii++)
		for (ij=0; ij<nr*NPEAK+NBACK; ij++)
			d[ii][ij] = 0.0;

	/* physicalize data position */
	w1 = inu*ks->delta_nu;
	tht = itht*TWOPI/ks->ntheta;
	akt = ks->k;

/* add background */
	back = p[offset]/(1.0+pow(w1/p[offset+1], p[offset+2]));
	backw = back*(1.0 + p[offset+3]*cos(2.0*(tht-p[offset+4])));
	back2 = 0.5*p[offset+5]*p[offset+7]
				/ ((w1-p[offset+6])*(w1-p[offset+6]) + 0.25*p[offset+7]*p[offset+7]);

	/* sum over ridges */
	for (ii=0; ii<nr; ii++)
	{
		A = p[ii*NPEAK+1];
		G = p[ii*NPEAK+2];
		twot = 2.0*(tht-p[ii*NPEAK+6]);
		co = cos(twot);
		anis = 1.0+p[ii*NPEAK+5]*co;
		shift = akt*(p[ii*NPEAK+3]*cos(tht) + p[ii*NPEAK+4]*sin(tht))/TWOPI;
		den1 = w1 - p[ii*NPEAK] + shift;
		den = den1*den1 + 0.25*G*G;
		lor = 0.5*A*G * anis / den;

		/* compute derivative and place in d[] */

		/* PRIMARY: 0 */
		d[ii*NPEAK+0][ii*NPEAK+0] = 4.*A*G*anis*den1*den1/(den*den*den) - A*G*anis/(den*den);
		d[ii*NPEAK+0][ii*NPEAK+1] = G*anis*den1/(den*den);
		d[ii*NPEAK+0][ii*NPEAK+2] = -2.*A*G*G*anis*den1/(den*den*den) + A*anis*den1/(den*den);
		d[ii*NPEAK+0][ii*NPEAK+3] = -2.*A*G*anis*den1*den1*akt*cos(tht)/(den*den*den*PI) + A*G*anis*akt*cos(tht)/(den*den*TWOPI);
		d[ii*NPEAK+0][ii*NPEAK+4] = -2.*A*G*anis*den1*den1*akt*sin(tht)/(den*den*den*PI) + A*G*anis*akt*sin(tht)/(den*den*TWOPI);
		d[ii*NPEAK+0][ii*NPEAK+5] = A*G*co*den1/(den*den);
		d[ii*NPEAK+0][ii*NPEAK+6] = 2.*A*G*p[ii*NPEAK+5]*sin(twot)*den1/(den*den);
		/* PRIMARY: 1 */
		d[ii*NPEAK+1][ii*NPEAK+1] = 0.0;
		d[ii*NPEAK+1][ii*NPEAK+2] = 0.5*anis/den - 0.5*G*G*anis/(den*den);
		d[ii*NPEAK+1][ii*NPEAK+3] = -G*anis*den1*akt*cos(tht)/(den*den*TWOPI);
		d[ii*NPEAK+1][ii*NPEAK+4] = -G*anis*den1*akt*sin(tht)/(den*den*TWOPI);
		d[ii*NPEAK+1][ii*NPEAK+5] = 0.5*G*co/den;
		d[ii*NPEAK+1][ii*NPEAK+6] = G*p[ii*NPEAK+5]*sin(twot)/den;
		/* PRIMARY: 2 */
		d[ii*NPEAK+2][ii*NPEAK+2] = -1.5*A*G*anis/(den*den) + A*G*G*G*anis/(den*den*den);
		d[ii*NPEAK+2][ii*NPEAK+3] = -A*anis*den1*akt*cos(tht)/(TWOPI*den*den) + 2.*A*G*G*anis*den1*akt*cos(tht)/(TWOPI*den*den*den);
		d[ii*NPEAK+2][ii*NPEAK+4] = -A*anis*den1*akt*sin(tht)/(TWOPI*den*den) + 2.*A*G*G*anis*den1*akt*sin(tht)/(TWOPI*den*den*den);
		d[ii*NPEAK+2][ii*NPEAK+5] = 0.5*A*co/den - 0.5*A*G*G*co/(den*den);
		d[ii*NPEAK+2][ii*NPEAK+6] = A*p[ii*NPEAK+5]*sin(twot)/den - A*G*G*p[ii*NPEAK+5]*sin(twot)/(den*den);
		/* PRIMARY: 3 */
		d[ii*NPEAK+3][ii*NPEAK+3] = -akt*akt*A*G*anis*cos(tht)*cos(tht)/(TWOPI*TWOPI*den*den) + 4.*akt*akt*A*G*anis*cos(tht)*cos(tht)*den1*den1/(TWOPI*TWOPI*den*den*den);
		d[ii*NPEAK+3][ii*NPEAK+4] = -akt*akt*A*G*anis*cos(tht)*sin(tht)/(TWOPI*TWOPI*den*den) + 4.*akt*akt*A*G*anis*cos(tht)*sin(tht)*den1*den1/(TWOPI*TWOPI*den*den*den);
		d[ii*NPEAK+3][ii*NPEAK+5] = -akt*A*G*co*cos(tht)*den1/(TWOPI*den*den);
		d[ii*NPEAK+3][ii*NPEAK+6] = -akt*A*G*p[ii*NPEAK+5]*sin(twot)*cos(tht)*den1/(PI*den*den);
		/* PRIMARY: 4 */
		d[ii*NPEAK+4][ii*NPEAK+4] = -akt*akt*A*G*anis*sin(tht)*sin(tht)/(TWOPI*TWOPI*den*den) + 4.*akt*akt*A*G*anis*sin(tht)*sin(tht)*den1*den1/(TWOPI*TWOPI*den*den*den);
		d[ii*NPEAK+4][ii*NPEAK+5] = -akt*A*G*co*sin(tht)*den1/(TWOPI*den*den);
		d[ii*NPEAK+4][ii*NPEAK+6] = -akt*A*G*p[ii*NPEAK+5]*sin(twot)*sin(tht)*den1/(PI*den*den);
		/* PRIMARY: 5 */
		d[ii*NPEAK+5][ii*NPEAK+5] = 0.0;
		d[ii*NPEAK+5][ii*NPEAK+6] = -A*G*sin(twot)/den;
		/*  PRIMARY: 6 */
		d[ii*NPEAK+6][ii*NPEAK+6] = -2.*A*G*co/den;
	}

	/* compute background */
	twot = 2.*(tht-p[offset+4]);
	co = cos(twot);
	anis = 1.0 + p[offset+3]*co;
	den1 = w1/p[offset+1];
	denp = pow(den1, p[offset+2]);
	den = denp+1.;
	po = p[offset]*anis/den;
	/* PRIMARY: 0 */
	d[offset+0][offset+0] = 0.0;
	d[offset+0][offset+1] = p[offset+2]*denp*po/den/p[offset+1];
	d[offset+0][offset+2] = -denp*log(den)*anis/(den*den);
	d[offset+0][offset+3] = co/den;
	d[offset+0][offset+4] = 2.*p[offset+3]*sin(twot)/den;
	/* PRIMARY: 1 */
	d[offset+1][offset+1] = -(p[offset+2]+1.)*denp*po/(p[offset+1]*p[offset+1]*den) + 2.*p[offset+2]*p[offset+2]*denp*denp*po/(den*den*p[offset+1]*p[offset+1]);
	d[offset+1][offset+2] = denp*po/(den*p[offset+1]) + p[offset+2]*denp*log(den1)*po/(p[offset+1]*den) - 2.*p[offset+2]*denp*denp*log(den1)*po/(p[offset+1]*den*den);
	d[offset+1][offset+3] = p[offset+2]*denp*p[offset]*anis/(p[offset+1]*den*den);
	d[offset+1][offset+4] = 2.*p[offset+2]*denp*p[offset]*p[offset+3]*sin(twot)/(p[offset+1]*den*den);
	/* PRIMARY: 2 */
	d[offset+2][offset+2] = -denp*log(den1)*log(den1)*po/den + 2.*denp*denp*log(den1)*log(den1)*po/(den*den);
	d[offset+2][offset+3] = -denp*log(den1)*p[offset]*co/(den*den);
	d[offset+2][offset+4] = -2.*denp*log(den1)*p[offset]*p[offset+3]*sin(twot)/(den*den);
	/* PRIMARY: 3 */
	d[offset+3][offset+3] = 0.0;
	d[offset+3][offset+4] = 2.*p[offset]*sin(twot)/den;
	/* PRIMARY: 4 */
	d[offset+4][offset+4] = -4.*p[offset]*(anis-1.)/den;

	/* background lorentzian */
	A = p[offset+5];
	G = p[offset+7];
	den1 = w1 - p[offset+6];
	den = den1*den1 + 0.25*G*G;
	lor = 0.5*A*G / den;
	/* PRIMARY: 5 amp*/
	d[offset+5][offset+5] = 0.0;
	d[offset+5][offset+6] = G*den1/(den*den);
	d[offset+5][offset+7] = 0.5/den - 0.5*G*G/(den*den);
	/* PRIMARY: 6 freq */
	d[offset+6][offset+6] = -A*G/(den*den) - 4.*A*G*den1*den1/(den*den*den);
	d[offset+6][offset+7] = A*den1/(den*den) - 2.*A*G*G*den1/(den*den*den);
	/* PRIMARY: 7 width */
	d[offset+7][offset+7] = -1.5*A*G/(den*den) + A*G*G*G/(den*den*den);
	
	/* Ensure symmetry */
	for (ii=0; ii<nr*NPEAK+NBACK; ii++)
		for (ij=ii; ij<nr*NPEAK+NBACK; ij++)
			d[ij][ii] = d[ii][ij];
}

/* compute all first derivatives at a single point */
void dm (double *p, struct kslice *ks, int inu, int itht, double *d)
{
	int ii, ij, nr, offset;
	double lor, w1, tht, akt, back, backw, back2, den, co, A, G, twot, anis, shift, den1;
	double denp, po;

	nr = ks->n;
	offset = nr*NPEAK;

	/* physicalize data position */
	w1 = inu*ks->delta_nu;
	tht = itht*TWOPI/ks->ntheta;
	akt = ks->k;

/* add background */
	back = p[offset]/(1.0+pow(w1/p[offset+1], p[offset+2]));
	backw = back*(1.0 + p[offset+3]*cos(2.0*(tht-p[offset+4])));
	back2 = 0.5*p[offset+5]*p[offset+7]
				/ ((w1-p[offset+6])*(w1-p[offset+6]) + 0.25*p[offset+7]*p[offset+7]);

	/* loop over ridges */
	for (ii=0; ii<nr; ii++)
	{
		A = p[ii*NPEAK+1];
		G = p[ii*NPEAK+2];
		twot = 2.0*(tht-p[ii*NPEAK+6]);
		co = cos(twot);
		anis = 1.0+p[ii*NPEAK+5]*co;
		shift = akt*(p[ii*NPEAK+3]*cos(tht) + p[ii*NPEAK+4]*sin(tht))/TWOPI;
		den1 = w1 - p[ii*NPEAK] + shift;
		den = den1*den1 + 0.25*G*G;
		lor = 0.5*A*G * anis / den;

		/* compute derivative and place in d[] */
		d[ii*NPEAK+0] = A*G*anis*den1/(den*den);
		d[ii*NPEAK+1] = 0.5*G*anis/den;
		d[ii*NPEAK+2] = 0.5*A*anis/den - 0.5*A*G*G*anis/(den*den);
		d[ii*NPEAK+3] = -akt*A*G*anis*cos(tht)*den1/(TWOPI*den*den); 
		d[ii*NPEAK+4] = -akt*A*G*anis*sin(tht)*den1/(TWOPI*den*den);
		d[ii*NPEAK+5] = 0.5*A*G*co/den;
		d[ii*NPEAK+6] = A*G*p[ii*NPEAK+5]*sin(twot)/den; 
	}
	/* compute background */
	if (p[offset+0] > 0.0)
	{
		twot = 2.*(tht-p[offset+4]);
		co = cos(twot);
		anis = 1.0 + p[offset+3]*co;
		den1 = w1/p[offset+1];
		denp = pow(den1, p[offset+2]);
		den = denp+1.;
		po = p[offset]*anis/den;

		d[offset+0] = anis/den;
		d[offset+1] = p[offset+2]*denp*A*anis/(den*den*p[offset+1]); 
		d[offset+2] = -denp*log(den1)*po/den; 
		d[offset+3] = A*co/den;
		d[offset+4] = 2.*A*p[offset+3]*sin(twot)/den;
	} else {
		d[offset+0] = 0.0;
		d[offset+1] = 0.0;
		d[offset+2] = 0.0;
		d[offset+3] = 0.0;
		d[offset+4] = 0.0;

	}

	if (p[offset+5]>0.0)
	{
		A = p[offset+5];
		G = p[offset+7];
		den1 = w1 - p[offset+6];
		den = den1*den1 + 0.25*G*G;
		lor = 0.5*A*G / den;

		d[offset+5] = 0.5*G/den;
		d[offset+6] = A*G*den1/(den*den);
		d[offset+7] = 0.5*A/den - 0.5*A*G*G/(den*den);
	} else {
		d[offset+5] = 0.0;
		d[offset+6] = 0.0;
		d[offset+7] = 0.0;
	}
}

/* compute model at a single point */
double model (double *p, struct kslice *ks, int inu, int itht)
{
	int ii, ij, nr, offset;
	double m, w1, tht, akt, back, den;

	nr = ks->n;
	offset = nr*NPEAK;

	m = 0;

	/* physicalize data position */
	w1 = inu*ks->delta_nu;
	tht = itht*TWOPI/ks->ntheta;
	akt = ks->k;

	/* add background */
	back = p[offset]/(1.0+pow(w1/p[offset+1], p[offset+2]));
	back *= (1.0 + p[offset+3]*cos(2.0*(tht-p[offset+4])));
	back += 0.5*p[offset+5]*p[offset+7]
				/ ((w1-p[offset+6])*(w1-p[offset+6]) + 0.25*p[offset+7]*p[offset+7]);

	/* sum over ridges */
	for (ii=0; ii<nr; ii++)
	{
		den = w1 - p[ii*NPEAK] + akt*(p[ii*NPEAK+3]*cos(tht) + p[ii*NPEAK+4]*sin(tht))/TWOPI;
		den = den*den + 0.25*p[ii*NPEAK+2]*p[ii*NPEAK+2];
		m += 0.5*p[ii*NPEAK+1]*p[ii*NPEAK+2]*(1.+p[ii*NPEAK+5]*cos(2.*(tht-p[ii*NPEAK+6])))/den;
	}

	return back + m;
}

