#include "fitsio.h"
#include "mpfit.h"

#define NPEAK (7) /* number of parameters for each peak */
#define NBACK (6) /* number of parameters for background terms */

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

/* Set which values correspond to which modes */
#define WEIGHT_NOISE (0)
#define WEIGHT_MAXL (1)
#define NOISE_CONST (0)
#define NOISE_SMOOTH (1)
#define NOISE_WAVELET (2)


struct kslice
{
	float **data, **noise;
	int start, end, ntheta;
	float delta_nu, k;
	struct params* par;
};

struct params
{
	char *fitsfname, *modelfname, *outfname;
	int silent, chiweight, noisemode;
	int kstart, kend;
	double ftol, xtol, gtol;
	int niter;
	char *debugfname, *covarfname;
};

/* main.c */
double model (int numridges, double nu, double k, double theta, double* p);
int funk(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* noise.c */
void compute_noise_smooth (float*** pol, float*** noise, int nnu, int nk, int ntheta, int radius);
void compute_noise_const (float***pol, float*** noise, int nnu, int nk, int ntheta);
void compute_noise_wavelet (float*** pol, float*** noise, int nnu, int nk, int ntheta);
unsigned int powerof2 (unsigned int n);

/* fit.c */
int fit_peak (float *freq, float *amp, float *width, float*** pol, float*** noise, float delta_nu, float delta_k, int nnu, int ntheta, int k);
int funk_single(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int fit_back (double* amp, double* cutoff, double* power, float*** pol, float*** noise, float delta_nu, int nnu, int ntheta, int k);
int funk_back(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* io.c */
void output_debug (struct params* p);
void output_covar (double* covar, int n, struct params* p);
int read_fits_file (float**** spec, float**** noise, struct params* p, 
		int* ntheta, int* nk, int* nnu, float* delta_k, float* delta_nu);
void read_param_file (char* fname, struct params* p);
void trim (char* str);
