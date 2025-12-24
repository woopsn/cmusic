# MUSIC / MVDR Beamforming in C

This repository contains a C implementation of the MUSIC (Multiple Signal Classification) algorithm combined with MVDR (Minimum Variance Distortionless Response) beamforming.

MUSIC is a technique for estimating the direction of arrival (DOA) of multiple signals received at a sensor array. It can provide high-resolution DOA estimates for uncorrelated narrowband sources in the presence of noise.

**Dependencies**

  - -llapack and -lblas for linear algebra operations

## Overview - sensor array processing

The basic assumption behind MUSIC is that narrowband signals originating in the sensors' "far field" arrive at the array with a planar wavefront. The time delay between the wave's arrival at each sensor can then be represented as a phase shift `a_i(θ)` determined by the wave's direction of arrival `θ` the location of the sensor. The instantaneous snapshot of a signal with baseband representation `s(t)` is the vector `x = (a_1(θ), a_2(θ), ..., a_M(θ))^T s(t) + n(t)` - the signal lags in phase by some amount at each sensor and is combined with a noise process `n(t)`.

The term `a(θ) = (a_1(θ), a_2(θ), ..., a_M(θ))^T` is known as the "steering" vector, and as `θ` varies it traces a smooth curve in `ℂ^M` (the array manifold).

MUSIC is in the family of "subspace methods" for DOA estimation. `N` snapshots are collected and used to compute a covariance estimate `R = (1/N) X X^H`, which is then decomposed into an eigenbasis (using LAPACK/BLAS). Assuming there are `K < M` uncorrelated sources and uncorrelated noise, the matrix `R` will be full-rank and normal leading to `M` orthogonal eigenvectors. `M-K` of these correspond to noise only, and when the SNR is sufficiently high these subspaces are associated with the `M-K` smallest eigenvalues of `R`.

The "signal" and "noise" subspaces of `R` are thus defined to be the span of the `K` largest and `M-K` smallest eigenvectors respectively.

## Overview - the pseudospectrum

When SNR is sufficiently high, the signal subspace of `R` is approximately the span of steering vectors corresponding to actual DOAs of the `K` (uncorrelated) sources. These vectors lie on the array manifold, and the smoothness of the steering vector function implies that nearby `a(θ)` are highly correlated with the signal subspace as well. The projections of `a(θ)` onto the signal and noise subspaces varies smoothly with `θ`, resulting in local extrema near the true DOAs. Typically a "pseudospectrum" will be computed over some range of angles and then its peaks located numerically, although analytic/root finding methods may be applicable in some cases.

In MUSIC, the pseudospectrum at `θ` is the reciprocal of the projection of `a(θ)` onto the noise subspace, and the estimated DOAs are associated with the `K` largest peaks of this function. . . It's just that simple!

## Usage

```sh
# Simulate 2 independent sources, ULA with 8 sensors, 200 snapshots
$ ./music -K 2 -M 8 -N 200

# Estimates, pseudospectrum, and data matrix are logged as CSV
$ ls -l output
DOA.csv S_MUSIC.csv X.csv
```

## TODO

- Allow complication of the ULA and noise models
- MVDR beamforming implementation

