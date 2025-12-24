#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <stdio.h>

#include "ula.h"

double nrand(double mean, double stddev);

// Generate ULA sensor data for K sources given d.o.a. (T) and power (P) parameters
//
// TODO
//   - allow configurable sensor spacing d
//   - allow non-uniform noise power across sensors
//   - allow non-Gaussian signals/noise
ULAData* generate_ula_data(int M, int N, int K, double Pn, double* T, double *P) {
  ULAData *data = ula_data_alloc(M, 1/4.0, N); // assuming d=1/4 wavelengths

  // For each source, baseband signal arriving from T(heta)[k] with power P[k]
  // is rotated by steering vector and summed at each sensor.
  for (int k = 0; k < K; k++) {
    double theta = T[k] * M_PI / 180.0;

    for (int n = 0; n < N; n++) {
      // Snapshot baseband signal
      double s_re = nrand(0.0, sqrt(P[k]/2.0));
      double s_im = nrand(0.0, sqrt(P[k]/2.0));

      for (int m = 0; m < M; m++) {
        double x_mn_re, x_mn_im;

        // Apply phase offset and steering vector for sensor m
        double phase = 2.0 * M_PI * m * data->d * sin(theta);
        x_mn_re = s_re * cos(phase) - s_im * sin(phase);
        x_mn_im = s_re * sin(phase) + s_im * cos(phase);

        data->data[2*(m*N + n)] += x_mn_re;
        data->data[2*(m*N + n) + 1] += x_mn_im;
      }
    }
  }

  // Add wgn to ula with uniform power Pn
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      double noise_real = nrand(0.0, sqrt(Pn/2.0));
      double noise_imag = nrand(0.0, sqrt(Pn/2.0));
  
      data->data[2*(m*N + n)] += noise_real;
      data->data[2*(m*N + n) + 1] += noise_imag;
    }
  }

  return data;
}

double nrand(double mean, double stddev) {
  // Box-Muller transform
  double u1 = ((double)rand() / RAND_MAX);
  double u2 = ((double)rand() / RAND_MAX);
  double z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
  return z0 * stddev + mean;
}
