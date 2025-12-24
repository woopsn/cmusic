#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <lapacke.h>

#include "ula.h"
#include "music.h"

#define ANGLE_GRID_STEP 0.1

typedef lapack_complex_double pair;

pair* covariance_alloc(ULAData *data);
int svd_alloc(pair **w, pair **vr, int M, pair *R);
MusicResult* spectrum_alloc(int M, pair *w, pair *vr, int K);
void spectrum_analyze(MusicResult *res, int K);

/*
 * 1. Compute covariance matrix R
 * 2. Compute eigendecomposition of R
 * 3. Separate K-dimensional signal subspace from noise subspace assuming SNR >> 1
 * 4. For each angle \theta, compute the MUSIC spectrum
 * 5. DOA estimates = argument of top K peaks
 */
MusicResult* music(ULAData *data, int K) {
  int M = data->M; // number of sensors
  pair *R = covariance_alloc(data);

  // Compute eigendecomposition of R using LAPACK
  pair *w, *vr;
  svd_alloc(&w, &vr, M, R);

  // Calculate MUSIC spectrum
  MusicResult *res = spectrum_alloc(M, w, vr, K);

  // Find largest K peaks
  spectrum_analyze(res, K);

  return res;
}

void music_result_free(MusicResult *res) {
  free(res->spectrum);
  free(res->angles);
  free(res);
}

pair* covariance_alloc(ULAData *data) {
  int M = data->M;
  int N = data->N;

  // Allocate matrix and initialize to 0
  pair *R = (pair *)malloc(M * M * sizeof(pair));
  for (int i = 0; i < M * M; i++)
    R[i] = 0.0 + I*0.0;

  // Compute sample covariance R = (1/N) * X * X^H
  for (int i = 0; i < M; i++) {
    for (int n = 0; n < N; n++) {
      pair xi = data->data[2*(i*N + n)] + I*data->data[2*(i*N + n) + 1];
      for (int j = 0; j < M; j++) {
        pair xj = data->data[2*(j*N + n)] + I*data->data[2*(j*N + n) + 1];
        R[i*M + j] += xi * conj(xj);
      }
    }
  }

  double norm_factor = 1.0 / N;
  for (int i = 0; i < M * M; i++)
    R[i] *= norm_factor;

  return R;
}

int svd_alloc(pair **w, pair **vr, int M, pair *R) {
  lapack_int n = M;
  lapack_int lda = M;

  pair *vl = NULL; // left eigenvectors not needed
  pair *eigs = (pair *)malloc(M * sizeof(pair)); // eigenvalues
  pair *eigvecs = (pair *)malloc(M * M * sizeof(pair)); // right eigenvectors

  lapack_int info = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, R, lda, eigs, vl, lda, eigvecs, lda);

  if (info != 0) {
    fprintf(stderr, "Error in eigen decomposition: %d\n", info);
    free(eigs);
    free(eigvecs);
    return -1;
  }

  // Simple bubble sort for eigenvalues and corresponding eigenvectors
  for (int i = 0; i < M - 1; i++) {
    for (int j = 0; j < M - i - 1; j++) {
      if (creal(eigs[j]) < creal(eigs[j + 1])) {
        // Swap eigenvalues
        pair temp_eig = eigs[j];
        eigs[j] = eigs[j + 1];
        eigs[j + 1] = temp_eig;

        // Swap corresponding eigenvectors
        for (int k = 0; k < M; k++) {
          pair temp_eigvec = eigvecs[k*M + j];
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

MusicResult* spectrum_alloc(int M, pair *w, pair *vr, int K) {
  int n_angles = (int)(180.0 / ANGLE_GRID_STEP) + 1;

  MusicResult *res = (MusicResult *)malloc(sizeof(MusicResult));
  res->n_angles = n_angles;
  res->spectrum =  (double *)malloc(n_angles * sizeof(double));
  res->n_peaks = K;
  res->angles = (double *)malloc(K * sizeof(double));

  // Calculate the pseudo spectrum
  // S_MUSIC(theta) = 
  pair a_theta[M]; // Steering vector a(theta)
  for (int i = 0; i < n_angles; i++) {
    double theta = -90.0 + i * ANGLE_GRID_STEP;
    double theta_rad = theta * M_PI / 180.0;

    for (int m = 0; m < M; m++) {
      double phase = 2.0 * M_PI * m * 0.25 * sin(theta_rad); // assuming d=1/4 wavelengths
      a_theta[m] = cos(phase) + I * sin(phase);
    }

    // Project onto noise subspace
    pair denominator = 0.0 + I*0.0;

    // For each noise vector (indices K to M-1)
    // compute |v_noise^H * a(theta)|^2 and add to denominator
    for (int k = K; k < M; k++) {
      pair v_noise_dot_a = 0.0 + I*0.0;
      for (int m = 0; m < M; m++)
        v_noise_dot_a += conj(vr[m*M + k]) * a_theta[m];

      denominator += v_noise_dot_a * conj(v_noise_dot_a);
    }

    res->spectrum[i] = 1.0 / creal(denominator);
  }

  return res;
}

void spectrum_analyze(MusicResult *res, int K) {
  int n_angles = res->n_angles;

  for (int k = 0; k < K; k++) {
    double max_value = -1.0;
    int max_index = -1;

    for (int i = 1; i < n_angles - 1; i++) {
      if (res->spectrum[i] > res->spectrum[i - 1] &&
          res->spectrum[i] > res->spectrum[i + 1] &&
          res->spectrum[i] > max_value) {
        // exclude previously found peaks
        for (int j = 0; j < k; j++) {
          if (fabs((-90.0 + i * ANGLE_GRID_STEP) - res->angles[j]) < ANGLE_GRID_STEP)
            goto skip_peak;
        }

        max_value = res->spectrum[i];
        max_index = i;
      }
    skip_peak:;
    }

    if (max_index != -1) {
      res->angles[k] = -90.0 + max_index * ANGLE_GRID_STEP;
    } else {
      res->angles[k] = 0.0; // default value if no peak found
    }
  }
}
