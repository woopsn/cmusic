# MUSIC / MVDR Beamforming in C

This repository contains a C implementation of the MUSIC (Multiple Signal Classification) algorithm combined with MVDR (Minimum Variance Distortionless Response) beamforming.

MUSIC is a technique for estimating the direction of arrival (DOA) of multiple signals received at a sensor array. It can provide high-resolution DOA estimates for uncorrelated narrowband sources in the presence of noise.

**Dependencies**

  - -llapack and -lblas for linear algebra operations

## Usage

The program `music` reads sensor array from an input CSV file, computes the MUSIC pseudospectrum, and estimates the DOAs of the sources. The results are written in CSV format to output/DOA.csv and output/S_MUSIC.csv.

The `-K`, `-M`, and `-N` flags are required in all cases.

    - `-K <int>`: Number of sources
    - `-M <int>`: Number of sensors
    - `-N <int>`: Number of snapshots

All angles are in degrees measured from the broadside of the array, with 0° corresponding to a signal arriving perpendicular to it and 90° to one that is completely aligned.

When input data is provided with `-i X.csv`, the dimensions of the table must match those specified by `-M` and `-N`. Complex entries should be formatted `a+bi`, with `a` and `b` floating point numbers.

``` sh
# Input data from X.csv, 2 sources, 200 snapshots of 8 sensors
$ ./music -i X.csv -K 2 -M 8 -N 200

# Estimates, pseudospectrum, and data matrix are logged as CSV
$ ls -l output
DOA.csv S_MUSIC.csv
```

If no input data is specified, synthetic data will be generated for testing according to the specified parameters.

```sh
# Simulate 2 sources with random DOA, ULA with 8 sensors, 200 snapshots
$ ./music -K 2 -M 8 -N 200

# Simulate 2 sources, one at 30°, one at random angle
$ ./music -K 2 -M 8 -N 200 S1=30 S2=70
```

## TODO

- Improve the ULA and noise models, add complications
- MVDR beamforming implementation
