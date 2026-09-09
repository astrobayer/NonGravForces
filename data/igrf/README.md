# data/igrf — IGRF-14 coefficients (for WP-23, LP-3b)

`igrf14coeffs.txt`: the 14th generation International Geomagnetic Reference Field,
Schmidt semi-normalised spherical-harmonic coefficients (nT), main field 1900.0
to 2025.0 and secular variation 2025-2030, as published by IAGA V-MOD via NOAA
NCEI: https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt, downloaded
2026-09-07 (owner's permission of the same day), 42 411 bytes,
SHA-256 8f8d88403028fc4ee92c4f38d97b46e0a87e2cfc496045b43c9e26c1d6b0903c.

Read by `NonGravForces.io.igrfCoefficients` since 2026-09-09 (WP-23); the file
fetched that day is **byte-identical** to the digest above, so upstream has not
moved. `setupNonGravForces` fetches it when it is missing.

**A correction to this file's earlier claim:** the reference Fortran `igrf14.f`
is NOT at the same site. `https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f` answers
`200` with an HTML error page, which is a soft 404 — a status code is not a
file. No published IGRF-14 test-value set was found; the IAGA `geomag70`
package carries sample input and output for IGRF-**13** only. So
`env.magneticField` will be verified against an independently written mirror,
as the force kernel is (WP-33), rather than against values that do not exist.
