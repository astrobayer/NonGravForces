# NonGravForces 0.8.0

Non-gravitational forces, torques and partials for LEO satellites (GRACE-FO reference case). Third concept, designed from scratch after the audits of `SatelliteNonGravToolbox` v1 and `NonGravToolbox` v0.1.1.

* Specification and plan: `HANDOVER.md` (single source of truth)
* Theory and implementation: `doc/NonGravForces_Theory_and_Implementation.pdf`
* Help browser: `setupNonGravForces`, then `doc NonGravForces`

```matlab
run('../referenceSystem/setupReferenceSystem.m'); run('../gracefoReading/setupGraceFOReading.m')
setupNonGravForces; runNgfTests
```
Status: 49 of 53 catalogued units delivered and verified in MATLAB, 4 documented stubs. **LP-7 is complete**; the four remaining stubs are blocked on data the project does not hold (a MAG1B file, a CERES granule, IGRF-14 coefficients). Next is LP-R, the release level.

Licensed MIT (see `LICENSE`). The vendored NRLMSISE-00 source under `mex/nrlmsise00/` is NOT covered by it and its redistribution terms are unsettled — read `THIRD-PARTY-NOTICES.md` before making this repository public. The authority is `doc/functionlist.csv`; `HANDOVER.md` Part A carries the level-by-level record.
