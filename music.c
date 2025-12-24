#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <lapacke.h>

#include "ula.h"
#include "music.h"

int music_sample_eigenstructure(lapack_complex_double **w, lapack_complex_double **vr, int M, lapack_complex_double *R);

int spectrum_analyze(double *DOA, int K, ULAData *data);

/*
 * 1. Compute eigendecomposition of data covariance matrix
 * 2. Separate K-dimensional signal subspace from noise subspace assuming SNR >> 1
 * 3. For each angle \theta, compute the MUSIC spectrum
 * 4. DOA estimates = argument of top K peaks
 */
int music(double *DOA, int K, ULAData *data) {
  // ensure covariance matrix is computed
  ula_sample_cov(data);

  lapack_complex_double *w, *vr;
  ula_sample_cov_eig(data, &w, &vr);
  ula_sample_smusic(data, w, vr, K);

  int k = spectrum_analyze(DOA, K, data);

  free(w);
  free(vr);

  return k;
}

void ula_sample_smusic(ULAData *data, lapack_complex_double *w, lapack_complex_double *vr, int K) {
  int M = data->M;
  lapack_complex_double a_theta[M]; // Steering vector a(theta)
  
  if (data->S_MUSIC == NULL) {
    data->n_angles = (int)180.0 / ANGLE_GRID_STEP + 1;
    data->S_MUSIC = (double *)malloc(data->n_angles * sizeof(double));
  }

  // Sweep theta from -90 to 90 degrees, calculating the magnitude of the steering vector's noise
  // component at each angle. S_MUSIC inverts this value, creating peaks at angles where the noise
  // component is minimal (ie where the steering vector aligns best with the signal subspace)
  int n_angles = data->n_angles;
  for (int i = 0; i < n_angles; i++) {
    double theta = (-90 + i*ANGLE_GRID_STEP) * M_PI/180.0;
    double phase_step = 2.0*M_PI*data->d*sin(theta);

    for (int m = 0; m < M; m++)
      a_theta[m] = cos(m*phase_step) + I*sin(m*phase_step);

    // The K weakest eigenvectors of R (vr indices K to M-1) form a basis for the noise subspace.
    // Compute |vnoise_k^H * a(theta)|^2 for each of them and collect, to get proj_noise(a(theta))
    double noise = 0.0 + I*0.0;
    for (int k = K; k < M; k++) {
      lapack_complex_double noise_k = 0.0 + I*0.0;
      for (int m = 0; m < M; m++)
        noise_k += conj(vr[m*M + k]) * a_theta[m];
      
      noise += creal(noise_k * conj(noise_k));
    }

    data->S_MUSIC[i] = 1.0 / noise;
  }
}

// DOA array is filled with angles in degrees, estimated from the K biggest peaks in S_MUSIC
int spectrum_analyze(double *DOA, int K, ULAData *data) {
  int n_angles = data->n_angles;

  for (int k = 0; k < K; k++) {
    double max_value = -1.0;
    int max_index = -1;

    for (int i = 1; i < n_angles - 1; i++) {
      double theta = -90.0 + i * ANGLE_GRID_STEP;

      if (data->S_MUSIC[i] > data->S_MUSIC[i - 1] &&
          data->S_MUSIC[i] > data->S_MUSIC[i + 1] &&
          data->S_MUSIC[i] > max_value) {
        // exclude previously found peaks
        for (int j = 0; j < k; j++)
          if (fabs(theta - DOA[j]) < ANGLE_GRID_STEP)
            goto skip_peak;

        max_value = data->S_MUSIC[i];
        max_index = i;
      }
    skip_peak:;
    }

    if (max_index != -1) {
      DOA[k] = -90.0 + max_index*ANGLE_GRID_STEP;
    } else {
      DOA[k] = 0.0; // default value if no peak found
    }
  }

  return K;
}
