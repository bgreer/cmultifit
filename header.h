#include "fitsio.h"
#include "mpfit.h"

/* For later, optimizing */
#define NNU (1152)
#define NK (192)
#define NTHETA (20)

#define NPEAK (7) /* number of parameters for each peak */
#define NBACK (8) /* number of parameters for background terms */

/* Set order of parameters in parameter file */
#define PARAM_SPECTRUM (0)
#define PARAM_MODEL (1)
#define PARAM_OUTPUT (2)
#define PARAM_SILENT (3)
#define PARAM_KRANGE (4)
#define PARAM_FTOL (5)
#define PARAM_XTOL (6)
#define PARAM_GTOL (7)
#define PARAM_NITER (8)
#define PARAM_DEBUG (9)
#define PARAM_COVAR (10)
#define PARAM_BACK (11)
#define PARAM_FIT (12)
#define PARAM_CUTOFF (13)
#define PARAM_FITABOVE (14)
#define PARAM_RIDGE (15)

/* Set which values correspond to which modes */
#define WEIGHT_NOISE (0)
#define WEIGHT_MAXL (1)
#define NOISE_CONST (0)
#define NOISE_SMOOTH (1)
#define NOISE_WAVELET (2)

#define PI    (3.14159265358979324)
#define TWOPI (6.28318530717958648)

#define gvg gsl_vector_get
#define map(x, y) (erf((x)/(y)))

struct kslice
{
	double **data;
	int start, end, ntheta, n;
	double delta_nu, k;
	struct params* par;
};

struct params
{
	char *fitsfname, *modelfname, *outfname;
	int silent, dofits, fit_above, detect_ridge;
	double ac_cutoff;
	int kstart, kend;
	double ftol, xtol, gtol;
	int niter;
	char *debugfname, *covarfname, *backfname;
};

double *thtarr, *thtpow;

/* main.c */
void mpreturn_translate (int mpreturn);

/* fit.c */
int fit_peak (struct params* p, double *freq, double *amp, double *width, 
	double*** pol, double delta_nu, double delta_k, int nnu, int ntheta, int k);
int funk_single(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int fit_back (struct params* p, double* amp, double* cutoff, double* power, 
	double*** pol, double delta_nu, int nnu, int ntheta, int k);
int funk_back(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* io.c */
void output_debug (struct params* p, double***pol, int ntheta, 
	int nk, int nnu, int k, int m, int n, double* x, double delta_nu, double delta_k);
void output_matrix (double* covar, int n, struct params* p);
int read_model_file (struct params *par, int nk, int **numridges, double ***freq, double ***amp, double ***width);
int read_fits_file (double**** spec, struct params* p, 
		int* ntheta, int* nk, int* nnu, double* delta_k, double* delta_nu);
void read_param_file (char* fname, struct params* p);
void trim (char* str);

/* function.c */
int funk(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

/* errbars.c */
void errbars (int numparams, double *p, struct kslice *ks, double **covar);
void fisher (int numparams, double *p, struct kslice *ks, double *covar);
double model (double *p, struct kslice *ks, int inu, int itht);
void dm (double *p, struct kslice *ks, int inu, int itht, double *d);
void d2m (double *p, struct kslice *ks, int inu, int itht, double **d);

