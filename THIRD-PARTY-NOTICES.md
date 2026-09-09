# Third-party notices

NonGravForces is MIT licensed (see `LICENSE`). That covers the code written for
this toolbox. It does not cover the third-party source vendored below, and this
file records what that source is and what is known about the terms it carries.

## mex/nrlmsise00/ — NRLMSISE-00, C release

**Files:** `nrlmsise-00.c`, `nrlmsise-00.h`, `nrlmsise-00_data.c`,
`nrlmsise-00_test.c`, `DOCUMENTATION`, `README.md`. Vendored **unmodified**;
`NonGravForces.buildMex` compiles them and nothing in this toolbox patches them.
A digest of every source file is written beside the built binary, so a silent
edit would turn the MEX stale on the next `buildMex('Check', true)`.

**Model authors:** Mike Picone, Alan Hedin and Doug Drob (Naval Research
Laboratory), building on MSISE90.
**C release maintainer:** Dominik Brodowski (mail@brodo.de).

**The redistribution terms are NOT settled, and this is the open point.** The
only legal text in the package is section 1 of its `DOCUMENTATION`, which is a
warranty disclaimer plus a request to send patches and bug fixes back to the
maintainer. It grants no explicit right to redistribute, and none is implied by
its absence. NRLMSISE-00 itself originates with a US government laboratory,
whose work is frequently in the public domain, but *frequently* is not a legal
finding and this file will not pretend otherwise.

**Before this repository is made public**, one of the following has to happen
(HANDOVER.md obligation O-34, and it is the owner's to settle):

1. obtain the maintainer's written permission to redistribute the C release, and
   record it here; or
2. establish that the source is public domain as a work of the US government,
   and record the basis for that; or
3. remove `mex/nrlmsise00/` from the repository and have
   `NonGravForces.buildMex` fetch or point at a user-supplied copy, which costs
   the reader a download and costs this repository nothing legally; or
4. replace the density backend. `env.density` already carries an `exponential`
   backend and refuses rather than degrades when `msis00` is unavailable, so the
   toolbox does not stop working without it — only the fidelity tier does.

Option 3 is the cheapest and does not depend on anyone answering an email.

## Data files under data/

`data/eop/`, `data/spaceweather/` and `data/samples/` hold small excerpts of
IERS, CelesTrak and NASA/JPL GRACE-FO Level-1B products, kept for testing. They
are redistributed under the terms their providers publish; they are not covered
by the MIT licence above and are not this toolbox's to relicense.
