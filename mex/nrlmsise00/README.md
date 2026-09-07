# mex/nrlmsise00 — vendored NRLMSISE-00 C source (unmodified)

The C implementation of NRLMSISE-00 (Picone, Hedin, Drob, Aikin 2002) by
Dominik Brodowski, release 20041227 with the 20070727 fix (see `DOCUMENTATION`
§6), taken over byte for byte on 2026-09-07 from the owner's copy in
`NonGravToolbox/external/nrlmsise00`. The model authors released the FORTRAN
original as a US Government work; the C port is distributed by its author
without warranty (`DOCUMENTATION` §1). Nothing in this folder is edited by
NonGravForces: the only code of ours is `../ngfMsisGateway.c`.

SHA-256 of the vendored files (recorded so that a silent edit is detectable):

| File | SHA-256 |
|---|---|
| `nrlmsise-00.c` | 09692a00305486c88a024c50031f7f5882fdfbe8582a51af8add7c37e12cd760 |
| `nrlmsise-00.h` | f567fd1e405b27f33018391d4333c91248777627dc42e60f6e7d58c5dbcb395a |
| `nrlmsise-00_data.c` | d0b3022f3c3e7ffdf703cc0e0d02339dfee01c0bf1520b29c68d5f4afd8784d5 |
| `nrlmsise-00_test.c` | dce37c3c13f07632673695fd6b762590bc66ce34248407e0882b325640c1ef16 |
| `DOCUMENTATION` | ddeb6c430007af087acceaa422d4a5b15bcbfa84258383f6fd75888ccb82f920 |

`nrlmsise-00_test.c` is not compiled into the MEX; it is kept because its
input block defines the 17 published test cases whose outputs (`DOCUMENTATION`
§7) are the oracle of `tests/tNgfEnv.m`. `NonGravForces.buildMex` compiles
`nrlmsise-00.c` and `nrlmsise-00_data.c` with the gateway into
`+NonGravForces/+env/private/msisGateway.<mexext>`; the binary is never
committed (`.gitignore`).
