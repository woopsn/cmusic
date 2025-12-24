#include <stdlib.h>

#include "ula.h"

ULAData* ula_alloc(int M, double d, int N) {
  ULAData *array = (ULAData *)malloc(sizeof(ULAData));
  array->M = M;
  array->d = d;

  array->N = N;
  array->X = (double *)calloc(2*M*N, sizeof(double));

  array->R = NULL;
  array->S_MUSIC = NULL;
  array->n_angles = 0;

  return array;
}

// Compute the sample covariance matrix R = (1/N) * X * X^H
void ula_sample_cov(ULAData *data) {
  int M = data->M;
  int N = data->N;
  lapack_complex_double *R = data->R;

  if (R != NULL) {
    for (int i = 0; i < M*M; i++)  // Zero out existing covariance matrix
      R[i] = 0.0 + I*0.0;
  } else {
    R = (lapack_complex_double *)calloc(M*M, sizeof(lapack_complex_double));
    data->R = R;
  }

  // Compute sample covariance R = (1/N) * X * X^H
  for (int i = 0; i < M; i++) {
    for (int n = 0; n < N; n++) {
      lapack_complex_double xi = data->X[2*(i*N + n)] + I*data->X[2*(i*N + n) + 1];
      for (int j = 0; j < M; j++) {
        lapack_complex_double xj = data->X[2*(j*N + n)] + I*data->X[2*(j*N + n) + 1];
        R[i*M + j] += xi * conj(xj);
      }
    }
  }

  double norm_factor = 1.0 / N;
  for (int i = 0; i < M*M; i++)
    R[i] *= norm_factor;
}

int ula_sample_cov_eig(ULAData *data, lapack_complex_double **w, lapack_complex_double **vr) {
  int M = data->M;

  lapack_complex_double *eigs = (lapack_complex_double *)malloc(M*sizeof(lapack_complex_double)); // eigenvalues
  lapack_complex_double *eigvecs = (lapack_complex_double *)malloc(M*M*sizeof(lapack_complex_double)); // right eigenvectors

  // Compute the right eigenstructure of data->R
  // https://www.netlib.org/lapack/explore-html/d4/d68/group__geev_ga8656466cd9da86a0b3b6966e49a7ce40.html
  lapack_int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', M, data->R,
                                  M, eigs, NULL, M, eigvecs, M);

  if (info != 0) {
    // fprintf(stderr, "Error in eigen decomposition: %d\n", info);
    free(eigs);
    free(eigvecs);
    return -1;
  }

  // Simple bubble sort for eigenvalues and corresponding eigenvectors
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < M - i - 1; j++) {
      if (creal(eigs[j]) < creal(eigs[j + 1])) {
        // Swap eigenvalues
        lapack_complex_double temp_eig = eigs[j];
        eigs[j] = eigs[j + 1];
        eigs[j + 1] = temp_eig;

        // Swap corresponding eigenvectors
        for (int k = 0; k < M; k++) {
          lapack_complex_double temp_eigvec = eigvecs[k*M + j];
          eigvecs[k * M + j] = eigvecs[k*M + j + 1];
          eigvecs[k * M + j + 1] = temp_eigvec;
        }
      }
    }
  }

  *w = eigs;
  *vr = eigvecs;

  return 0;
}

void ula_free(ULAData *data) {
  free(data->X);

  if (data->R != NULL)
    free(data->R);

  if (data->S_MUSIC != NULL)
    free(data->S_MUSIC);

  free(data);
}
