struct MusicResult {
  int n_angles;        // Number of angles in the pseudospectrum
  double *spectrum;    // MUSIC pseudospectrum (-90 to +90 degrees)

  int n_peaks;        // Number of detected peaks
  double *angles;     // Estimated DOA angles in degrees
};

typedef struct MusicResult MusicResult;

// Function to perform MUSIC DOA estimation
MusicResult* music(ULAData *data, int nsources);

// Function to free MusicResult structure
void music_result_free(MusicResult *res);
