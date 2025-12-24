#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <stdio.h>

#include "ula.h"

double nrand(double mean, double stddev);

ULAData* read_ula_data(char *path, int M, int N) {
  char buf[256];
  FILE *fp = fopen(path, "rb");
  if (fp == NULL) return NULL;

  ULAData *data = ula_alloc(M, 1/4.0, N); // assuming d=1/4 wavelengths

  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      float re, im;
      char op;

      if (n != 0 && (fscanf(fp, "%c", buf) != 1 || *buf != ',')) {
        fprintf(stderr, "Error reading comma from file (line %d, column %d)\n", m+1, n+1);
        goto parse_failure;
      }

      if (fscanf(fp, "%f", &re) != 1) {
        fprintf(stderr, "Error reading real part from file (line %d, column %d)\n", m+1, n+1);
        goto parse_failure;
      }

      if (fscanf(fp, " %c ", &op) != 1 || (op != '+' && op != '-')) {
        fprintf(stderr, "Error reading imaginary sign from file (line %d, column %d)\n", m+1, n+1);
        goto parse_failure;
      }

      if (fscanf(fp, "%fi", &im) != 1) {
        fprintf(stderr, "Error reading imaginary part from file (line %d, column %d)\n", m+1, n+1);
        goto parse_failure;
      }

      data->X[2*(m*N + n)] = (double)re;
      data->X[2*(m*N + n) + 1] = (double)(op == '-' ? -im : im);
    }

    if (fscanf(fp, "%1[\r\n]", buf) != 1 && m < M-1) {
      fprintf(stderr, "Error reading newline from file (line %d)\n", m);
      goto parse_failure;
    }
  }

  fclose(fp);
  return data;

  parse_failure:
  fscanf(fp, "%4s", buf);
  fprintf(stderr, "Next: '%s'\n", buf);
  ula_free(data);
  fclose(fp);
  return NULL;
}

// Generate ULA sensor data for K sources given d.o.a. (T) and power (P) parameters
//
// TODO
//   - allow configurable sensor spacing d
//   - allow non-uniform noise power across sensors
//   - allow non-Gaussian signals/noise
ULAData* generate_ula_data(int M, int N, int K, double Pn, double* T, double *P) {
  ULAData *data = ula_alloc(M, 1/4.0, N); // assuming d=1/4 wavelengths

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

        // Narrowband assumption => same signal is received at each sensor,
        // just phase-shifted according to angle of arrival
        double phase = 2.0 * M_PI * m * data->d * sin(theta);
        x_mn_re = s_re * cos(phase) - s_im * sin(phase);
        x_mn_im = s_re * sin(phase) + s_im * cos(phase);

        data->X[2*(m*N + n)] += x_mn_re;
        data->X[2*(m*N + n) + 1] += x_mn_im;
      }
    }
  }

  // Add wgn to ula with uniform power Pn
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      double noise_real = nrand(0.0, sqrt(Pn/2.0));
      double noise_imag = nrand(0.0, sqrt(Pn/2.0));
  
      data->X[2*(m*N + n)] += noise_real;
      data->X[2*(m*N + n) + 1] += noise_imag;
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
