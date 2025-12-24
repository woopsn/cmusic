#include <stdlib.h>

#include "ula.h"

ULAData* ula_data_alloc(int M, double d, int N) {
  ULAData *array = (ULAData *)malloc(sizeof(ULAData));
  array->M = M;
  array->d = d;

  array->N = N;
  array->data = (double *)calloc(2*M*N, sizeof(double));

  return array;
}

void ula_data_free(ULAData *data) {
  free(data->data);
  free(data);
}
