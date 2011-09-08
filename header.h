#include "fitsio.h"
#include "mpfit.h"

#define NPEAK (8) /* number of parameters for each peak */
#define NBACK (8) /* number of parameters for background terms */

/* Set order of parameters in parameter file */
#define PARAM_SPECTRUM (0)
#define PARAM_MODEL (1)
#define PARAM_OUTPUT (2)
#define PARAM_SILENT (3)
#define PARAM_WEIGHT (4)
#define PARAM_NOISE (5)
#define PARAM_KRANGE (6)
#define PARAM_FTOL (7)
#define PARAM_XTOL (8)
#define PARAM_GTOL (9)
#define PARAM_NITER (10)
#define PARAM_DEBUG (11)
#define PARAM_COVAR (12)
#define PARAM_BACK (13)

/* Set which values correspond to which modes */
#define WEIGHT_NOISE (0)
#define WEIGHT_MAXL (1)
#define NOISE_CONST (0)
#define NOISE_SMOOTH (1)
#define NOISE_WAVELET (2)

#define PI    (3.14159265358979324)
#define TWOPI (6.28318530717958648)

struct kslice
{
	double **data, **noise;
	int start, end, ntheta;
	double delta_nu, k;
	struct params* par;
};

struct params
{
	char *fitsfname, *modelfname, *outfname;
	int silent, chiweight, noisemode;
	int kstart, kend;
	double ftol, xtol, gtol;
	int niter;
	char *debugfname, *covarfname, *backfname;
};

double *thtarr, *thtpow, *rdeviates, *ideviates;

/* main.c */
double model (int numridges, double nu, double k, double theta, double* p);
void normalize (double ****spec, double** norm, int nnu, int nk, int ntheta);

/* noise.c */
void compute_noise_smooth (double*** pol, double*** noise, int nnu, int nk, int ntheta, int radius);
void compute_noise_const (double***pol, double*** noise, int nnu, int nk, int ntheta);
void compute_noise_wavelet (double*** pol, double*** noise, int nnu, int nk, int ntheta);
unsigned int powerof2 (unsigned int n);

/* fit.c */
int fit_peak (struct params* p, double *freq, double *amp, double *width, 
	double*** pol, double*** noise, double delta_nu, double delta_k, int nnu, int ntheta, int k);
int funk_single(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int fit_back (struct params* p, double* amp, double* cutoff, double* power, 
	double*** pol, double*** noise, double delta_nu, int nnu, int ntheta, int k);
int funk_back(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* io.c */
void output_debug (struct params* p, double***pol, double***noise, int ntheta, 
	int nk, int nnu, int k, int m, int n, double* x, double delta_nu, double delta_k);
void output_covar (double* covar, int n, struct params* p);
int read_fits_file (double**** spec, double**** noise, struct params* p, 
		int* ntheta, int* nk, int* nnu, double* delta_k, double* delta_nu);
void read_param_file (char* fname, struct params* p);
void trim (char* str);

/* function.c */
int funk(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int funk_corr(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
void calc_derivs (int m, int n, double* p, double *deviates, double **derivs, void *private_data);
