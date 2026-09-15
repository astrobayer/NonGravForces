# Third-party notices

NonGravForces is MIT licensed (see `LICENSE`). That covers the code written for
this toolbox. It does not cover the third-party source fetched below, and this
file records what that source is and what is known about the terms it carries.

## NRLMSISE-00 — fetched, not redistributed (settled 2026-09-09)

This repository does **not** contain the NRLMSISE-00 C source, and the reason is
that the package carries no right to redistribute it. Its only legal text is
section 1 of its own `DOCUMENTATION`: a warranty disclaimer plus a request to
send patches to the maintainer. No redistribution grant, and none implied by its
absence. NRLMSISE-00 originates with a US government laboratory, whose work is
often public domain — but *often* is not a legal finding, and this file will not
make one.

So the source is obtained from its maintainer, and **`setupNonGravForces` does
that for you** the first time it runs:

    setupNonGravForces          % clones into mex/nrlmsise00 at the pinned commit

or by hand, if you prefer:

    git clone https://git.linta.de/~brodo/nrlmsise-00.git mex/nrlmsise00

The installer is the right place for a network call and a build step is not:
`NonGravForces.buildMex` looks in `'MsisSource'`, then `NGF_NRLMSISE00`, then
`mex/nrlmsise00`, and **refuses** if it finds nothing rather than fetching. A
build that reaches the network without being asked is a build nobody can audit.
Setup prints one line whenever it clones, whatever `'Verbose'` says, because a
network action has to be visible; `'FetchMsis', "never"` leaves the network
alone entirely, and `NGF_NRLMSISE00` suppresses the fetch because you already
have a copy. CI needs no step of its own — `ciGate` calls setup before it
builds.

**Model authors:** Mike Picone, Alan Hedin and Doug Drob (Naval Research
Laboratory), building on MSISE90. **C release maintainer:** Dominik Brodowski.

### Which release, and why it changed

Until 2026-09-09 this repository vendored a 2007-era copy. Two things were found
when the removal was carried out, and both are recorded because they shaped the
decision:

1. **That copy could not be reproduced from upstream.** Its `nrlmsise-00.c`
   digest matches no commit in the maintainer's entire history — twelve commits,
   all checked. It came from a local tree and may carry an applied patch or
   different line endings. Nobody but its holder could have obtained those exact
   bytes, which is what made pinning to it untenable for a public repository.
2. **The published oracle did not change.** Section 7 of the `DOCUMENTATION` —
   the 17-case testing-output table `tests/tNgfEnv` parses — is byte-identical
   between the old copy and upstream. Only the maintainer's email, the sign-off,
   and two release-note lines differ. The move therefore does not disturb the
   reference values.

The owner's ruling is to **track upstream and take its bug fixes** — since the
old copy, upstream has fixed Argon densities, a memory corruption, the gtd7
output at 32.5 km, and compiler warnings. CI pins the commit
(`bc9a2feba4344e74201281e563332688a4d09cc3`) rather than following a branch, so
that someone else's push cannot change this toolbox's density numbers with
nothing in its diff to explain it. Moving the pin is a deliberate act with a run
behind it.

The build digest covers whatever source is actually compiled, so a different
release marks the binary stale and forces a rebuild. That is the mechanism that
makes tracking upstream safe rather than merely convenient.

## NRLMSIS 2.1 — fetched, not redistributed (settled 2026-09-15)

> This software incorporates the MSIS® empirical atmospheric model software
> designed and provided by NRL. Use is governed by the Open Source Academic
> Research License Agreement contained in the file nrlmsis2.1_license.txt.

That statement is required by section 4c of the agreement from any package that
incorporates the software, and it is repeated here, in `mex/ngfMsis21Shim.F90`,
and by `setupNonGravForces` whenever it fetches the package. (In NRL's release
tarball the file is named `nrlmsis2.1_license..txt`, with two dots.)

**What the agreement allows, in its own words:** use, reproduction and
modification "solely for research, academic, and non-profit purposes"; "Any
commercial use is prohibited." The MIT licence of this toolbox does not and
cannot widen that: **the `msis21` density backend is available to you only on
NRL's terms.** Every other backend (`exponential`, `msis00`) is unaffected.

What this repository does about it:

- **Nothing of NRL's is committed.** `setupNonGravForces` downloads NRL's
  release tarball (`nrlmsis2.1.tar.gz` of 2022-04-04 from
  `map.nrl.navy.mil`), refuses it unless its SHA-256 is
  `41e47b29f795d36a5cc252b2858aa2a384c4a7323ace3d48d3ea2f2b37a1a6a8`, and
  unpacks it into `mex/nrlmsis21` (ignored by git). `'FetchMsis', "never"`
  leaves the network alone; `NGF_NRLMSIS21` or `buildMex('Msis21Source', ...)`
  point at a copy you already hold.
- **No NRL file is modified.** The seven model files are compiled as shipped
  (gfortran, `-O2 -cpp -DDBLE`) and called through two `bind(C)` wrappers of
  this toolbox's own (`mex/ngfMsis21Shim.F90`) and a C gateway
  (`mex/ngfMsis21Gateway.c`).
- **The agreement travels with the binary.** `NonGravForces.buildMex` copies
  `nrlmsis2.1_license..txt` beside `msis21Gateway.<mexext>` together with the
  parameter file, and the build digest covers both (4c: "Each copy of the
  Software ... shall include a file containing this Agreement").

**Model authors:** John Emmert, Douglas Drob and McArthur Jones Jr. (Naval
Research Laboratory). Reference: Emmert et al. (2021), NRLMSIS 2.0, *Earth and
Space Science* 8, e2020EA001321.

## Data files under data/

`data/eop/`, `data/spaceweather/`, `data/samples/` and `data/ceres/` hold small
excerpts of IERS, CelesTrak, NASA/JPL GRACE-FO Level-1B and NASA LaRC CERES EBAF
products, kept so the tests and demos have a real arc to run on. The CERES
excerpt is licensed by NASA LaRC under CC BY 4.0; its attribution and how it was
cut are in `data/ceres/README.md`. They are redistributed under the terms their
providers publish; they are not covered by the MIT licence above and are not
this toolbox's to relicense.

The toolbox is published **via git only** — `.mltbx` packaging was deprecated on
2026-09-09 — so the question of what a package may carry does not arise.
`setupNonGravForces` downloads these, only when they are missing, and otherwise
supplying them is the user's responsibility.

| source | fetched by setup | where from |
|---|---|---|
| Earth orientation | yes | `https://datacenter.iers.org/data/9/finals2000A.all` |
| space weather | yes, with one retry | `https://celestrak.org/SpaceData/SW-All.csv` |
| IGRF-14 coefficients | yes | `https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt` |
| CERES EBAF 4.2.1 (rsut, rlut, rsdt) | yes, 65 MB, SHA-256 checked | `https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon/` |
| GRACE-FO Level-1B | yes, 147 MB, members only | `https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/` |

Every endpoint was fetched before it was written down. IERS answered 200 with the
fixed-format file; CelesTrak answered 503 once and 200 on the retry, which is why
that fetch retries; the GFZ archive for 2020-01-01 holds all 65 products under
exactly the filenames this toolbox's fixtures use, so setup keeps only the three it
needs (GNV1B, SCA1B, MAG1B) and discards the rest. Since 2026-09-15 every fetched
file must also pass a check of its first lines before it replaces anything.

`'FetchData', "never"` leaves the network alone entirely.
