#include <lapacke.h>

#define ANGLE_GRID_STEP 0.1
#define N_ANGLES ((int)(180.0 / ANGLE_GRID_STEP) + 1)

// Struct representing a Uniform Linear Array (ULA) of sensors
struct ULAData {
  int M;              // Number of sensors
  double d;           // Spacing in wavelengths

  int N;              // Number of snapshots
  double *X;          // Data matrix - row major array, real and imaginary parts

  lapack_complex_double *R; // Data covariance matrix
  lapack_complex_double **w;  // Eigenvalues of R
  lapack_complex_double **vr; // Right eigenvectors of R

  int n_angles;             // Number of angles in the MUSIC spectrum
  double *S_MUSIC;          // MUSIC pseudospectrum
};

typedef struct ULAData ULAData;
typedef lapack_complex_double lapack_c;

#define ULA_X_REAL(data, m, n) data->X[2*((m)*(data->N) + (n))]
#define ULA_X_IMAG(data, m, n) data->X[2*((m)*(data->N) + (n)) + 1]

ULAData* ula_alloc(int M, double d, int N);

void     ula_sample_cov(ULAData *data);
int      ula_sample_cov_eig(ULAData *data, lapack_complex_double **w, lapack_complex_double **vr);

void     ula_sample_smusic(ULAData *data, lapack_complex_double *w, lapack_complex_double *vr, int K);

void     ula_free(ULAData *data);
