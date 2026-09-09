# Third-party notices

NonGravForces is MIT licensed (see `LICENSE`). That covers the code written for
this toolbox. It does not cover the third-party source vendored below, and this
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

## Data files under data/

`data/eop/`, `data/spaceweather/` and `data/samples/` hold small excerpts of
IERS, CelesTrak and NASA/JPL GRACE-FO Level-1B products, kept for testing. They
are redistributed under the terms their providers publish; they are not covered
by the MIT licence above and are not this toolbox's to relicense.
