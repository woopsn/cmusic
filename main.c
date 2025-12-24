#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <getopt.h>
#include <time.h>

#include "ula.h"
#include "music.h"
#include "gendata.h"

void parse_parameters(int *M, int *K, int *N, double *Pn, int argc, char** argv);
void parse_sources(int argc, char** argv, int K, double* T, double *P);
void write_log(ULAData *data, MusicResult *res, int M, int N, int K, double* T, double *P, double Pn);

int main(int argc, char** argv) {
  int M = 100;         // number of sensors
  int N = 200;         // number of snapshots
  int K = 3;           // number of sources. assume for now all have power (or variance)=1
  double Pn = 0.05;    // noise power
  parse_parameters(&M, &K, &N, &Pn, argc, argv);

  // srand(time(NULL));
  sranddev();

  double T[K]; // source DOAs in degrees (T=Theta)
  double P[K]; // source powers
  parse_sources(argc, argv, K, T, P);

  ULAData *data;
  MusicResult *result;

  data = generate_ula_data(M, N, K, Pn, T, P);
  result = music(data, K);

  write_log(data, result, M, N, K, T, P, Pn);

  music_result_free(result);
  ula_data_free(data);

  return 0;
}

void parse_parameters(int *M, int *K, int *N, double *Pn, int argc, char** argv) {
  int opt;
  while ((opt = getopt(argc, argv, "hM:K:N:s:")) != -1) {
    switch (opt) {
    case 'h':
    default:
      printf("Usage: %s [options] [specs ...]\n", argv[0]);
      printf("  -M <num_sensors>      Number of sensors (default: 100)\n");
      printf("  -K <num_sources>      Number of sources (default: 3)\n");
      printf("  -N <num_snapshots>    Number of snapshots (default: 200)\n\n");

      printf("Test data configuration:\n");
      printf("  -s <noise_power>      Noise power (default: 0.05)\n\n");

      printf("  Examples:\n");
      printf("    M=%d, K=%d, N=%d, Pn=%.2f:\n", *M, *K, *N, *Pn); // default values
      printf("      %s\n\n", argv[0]);

      printf("    S1 at 30 degrees, power 2.0; S2 at -15 degrees, power 1.0:\n");
      printf("      %s -K 2 S1=2.0@30.0 S2@-15.0\n\n", argv[0]);

      printf("    S1 at 30 degrees, power 2.0; S2 at random angle, power 1.0:\n");
      printf("      %s -K 2 S1=2.0@30.0\n\n", argv[0]);

      exit(opt == 'h' ? 0 : 1);

    case 'M':
      *M = atoi(optarg);
      break;

    case 'K':
      *K = atoi(optarg);
      break;

    case 'N':
      *N = atoi(optarg);
      break;

    case 's':
      *Pn = atof(optarg);
      break;
    }
  }
}

void parse_sources(int argc, char** argv, int K, double* T, double *P) {
  for (int k = 0; k < K; k++) {
    P[k] = 1.0; // all sources have power 1.0
    T[k] = -60.0 + ((double)rand() / RAND_MAX) * 120.0;
  }

  for (int i = optind; i < argc; i++) {
    char *arg = argv[i];
    if (arg[0] == 'S') {
      int src_idx = atoi(&arg[1]) - 1; // S1 -> index 0
      char *eq_ptr = strchr(arg, '=');
      char *at_ptr = strchr(arg, '@');

      if (eq_ptr != NULL) {
        double angle = atof(eq_ptr + 1);
        if (src_idx >= 0 && src_idx < K) {
          T[src_idx] = angle;
        }
      }

      if (at_ptr != NULL) {
        double power = atof(at_ptr + 1);
        if (src_idx >= 0 && src_idx < K) {
          P[src_idx] = power;
        }
      }
    }
  }
}

void write_log(ULAData *data, MusicResult *res, int M, int N, int K, double* T, double *P, double Pn) {
  FILE *f = fopen("output/X.csv", "w");
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      fprintf(f, "%lf+%lfi", data->data[2*(m*N + n)], data->data[2*(m*N + n) + 1]);
      if (n < N - 1) fprintf(f, ",");
    }
    fprintf(f, "\n");
  }
  fclose(f);

  f = fopen("output/S_MUSIC.csv", "w");
  for (int i = 0; i < res->n_angles; i++) {
    fprintf(f, "%lf,%lf\n", -90.0 + i * (180.0 / (res->n_angles - 1)), res->spectrum[i]);
  }
  fclose(f);

  f = fopen("output/DOA.csv", "w");
  for (int k = 0; k < K; k++) {
    fprintf(f, "%.2f, %.2f\n", T[k], res->angles[k]);
  }
}
