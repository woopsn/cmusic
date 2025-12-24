// Struct representing a Uniform Linear Array (ULA) of sensors
struct ULAData {
  int M;              // Number of sensors
  double d;           // Spacing in wavelengths

  int N;              // Number of snapshots
  double *data;       // 2D array of complex data (real and imag parts)
};

typedef struct ULAData ULAData;

ULAData *ula_data_alloc(int M, double d, int N);
void ula_data_free(ULAData *data);
