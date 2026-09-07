# data/igrf — IGRF-14 coefficients (for WP-23, LP-3b)

`igrf14coeffs.txt`: the 14th generation International Geomagnetic Reference Field,
Schmidt semi-normalised spherical-harmonic coefficients (nT), main field 1900.0
to 2025.0 and secular variation 2025-2030, as published by IAGA V-MOD via NOAA
NCEI: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt, downloaded
2026-09-07 (owner's permission of the same day), 42 411 bytes,
SHA-256 8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c.

Not read by any unit yet: `io.igrfCoefficients` and `env.magneticField` are
LP-3b (WP-23). The reference Fortran `igrf14.f` (same site) will be the oracle
source of published test values when that work package is written.
