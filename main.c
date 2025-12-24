#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <getopt.h>
#include <time.h>

#include "ula.h"
#include "music.h"
#include "gendata.h"

void parse_parameters(char **infile, int *M, int *K, int *N, double *Pn, int argc, char** argv);
void parse_sources(int argc, char** argv, int K, double* T, double *P);
void write_log(ULAData *data, double *DOA, int K, double* T, double *P, double Pn);

int main(int argc, char** argv) {
  int M = 100;         // number of sensors
  int N = 200;         // number of snapshots
  int K = 3;           // number of sources. assume for now all have power (or variance)=1
  double Pn = 0.05;    // noise power
  char *input_file = NULL;
  parse_parameters(&input_file, &M, &K, &N, &Pn, argc, argv);

  // srand(time(NULL));
  sranddev();

  ULAData *data;

  double T[K], P[K]; // synthetic source DOAs (in degrees) and power
  memset(T, 0, sizeof(T));
  memset(P, 0, sizeof(P));

  if (input_file != NULL) {
    data = read_ula_data(input_file, M, N);
    if (data == NULL) {
      fprintf(stderr, "Error reading input file %s\n", input_file);
      return 1;
    }
  } else {
    parse_sources(argc, argv, K, T, P);
    data = generate_ula_data(M, N, K, Pn, T, P);
  }

  double DOA_EST[K];
  music(DOA_EST, K, data);

  write_log(data, DOA_EST, K, T, P, Pn);
  ula_free(data);

  return 0;
}

void parse_parameters(char **infile, int *M, int *K, int *N, double *Pn, int argc, char** argv) {
  int opt;
  while ((opt = getopt(argc, argv, "i:M:K:N:s:h")) != -1) {
    switch (opt) {
    case 'h':
    default:
      printf("Usage: %s [options] [specs ...]\n", argv[0]);
      printf("  -M <num_sensors>      Number of sensors (default: 100)\n");
      printf("  -N <num_snapshots>    Number of snapshots (default: 200)\n");
      printf("  -i <input_file>       Input CSV containing sensor snapshot matrix\n");
      printf("                        M rows, N columns in a[+-]bi format\n\n");

      printf("If -i is not provided, data will be generated from the following parameters:\n");
      printf("  -K <num_sources>      Number of sources (default: 3)\n");
      printf("  -s <noise_power>      Noise power (default: 0.05)\n\n");

      printf("  Examples:\n");
      printf("    M=%d, K=%d, N=%d, Pn=%.2f:\n", *M, *K, *N, *Pn); // default values
      printf("      %s\n\n", argv[0]);

      printf("    S1 at 30 degrees, power 2.0; S2 at -15 degrees, power 1.0:\n");
      printf("      %s -K 2 S1=30.0@2.0 S2=-15.0\n\n", argv[0]);

      printf("    S1 at 30 degrees, power 2.0; S2 at random angle, power 1.0:\n");
      printf("      %s -K 2 S1=30.0@2.0\n\n", argv[0]);

      exit(opt == 'h' ? 0 : 1);

    case 'M':
      *M = atoi(optarg);
      break;

    case 'N':
      *N = atoi(optarg);
      break;

    case 'i':
      *infile = optarg;
      break;

    case 'K':
      *K = atoi(optarg);
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

void write_log(ULAData *data, double *DOA, int K, double* T, double *P, double Pn) {
  int M = data->M;
  int N = data->N;

  FILE *f = fopen("output/X.csv", "w");
  for (int m = 0; m < M; m++) {
    for (int n = 0; n < N; n++) {
      fprintf(f, "%lf+%lfi", data->X[2*(m*N + n)], data->X[2*(m*N + n) + 1]);
      if (n < N - 1) fprintf(f, ",");
    }
    fprintf(f, "\n");
  }
  fclose(f);

  f = fopen("output/S_MUSIC.csv", "w");
  for (int i = 0; i < data->n_angles; i++) {
    fprintf(f, "%lf,%lf\n", -90.0 + i*(180.0/(data->n_angles - 1)), data->S_MUSIC[i]);
  }
  fclose(f);

  f = fopen("output/DOA.csv", "w");
  for (int k = 0; k < K; k++) {
    fprintf(f, "%.2f, %.2f\n", T[k], DOA[k]);
  }
}
