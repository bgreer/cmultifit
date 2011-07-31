#include "fitsio.h"
#include "mpfit.h"

/* Set order of parameters in parameter file */
#define PARAM_SPECTRUM (0)
#define PARAM_MODEL (1)
#define PARAM_OUTPUT (2)
#define PARAM_SILENT (3)
#define PARAM_WEIGHT (4)
#define PARAM_NOISE (5)

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
};

struct params
{
	char *fitsfname, *modelfname, *outfname;
	int silent, chiweight, noisemode;
};

double model (int numridges, double nu, double k, double theta, double* p);
int funk(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

void compute_noise (float*** pol, float*** noise, int nnu, int nk, int ntheta, int radius);
void compute_noise_wavelet (float*** pol, float*** noise, int nnu, int nk, int ntheta, int radius);
unsigned int powerof2 (unsigned int n);

int fit_peak (float *freq, float *amp, float *width, float*** pol, float*** noise, float delta_nu, float delta_k, int nnu, int ntheta, int k);
int funk_single(int m, int n, double* p, double* deviates, double**derivs, void* private_data);
int fit_back (double* amp, double* cutoff, double* power, float*** pol, float*** noise, float delta_nu, int nnu, int ntheta, int k);
int funk_back(int m, int n, double* p, double* deviates, double**derivs, void* private_data);

void read_param_file (char* fname, struct params* p);
void trim (char* str);
