# NonGravForces 0.6.0

Non-gravitational forces, torques and partials for LEO satellites (GRACE-FO reference case). Third concept, designed from scratch after the audits of `SatelliteNonGravToolbox` v1 and `NonGravToolbox` v0.1.1.

* Specification and plan: `HANDOVER.md` (single source of truth)
* Theory and implementation: `doc/NonGravForces_Theory_and_Implementation.pdf`
* Help browser: `setupNonGravForces`, then `doc NonGravForces`

```matlab
run('../referenceSystem/setupReferenceSystem.m'); run('../gracefoReading/setupGraceFOReading.m')
setupNonGravForces; runNgfTests
```
Status: 43 of 52 units delivered and verified in MATLAB (core, data I/O, environment, geometry and attitude, kernel, quality control), 12 documented stubs. The authority is `doc/functionlist.csv`; `HANDOVER.md` Part A carries the level-by-level record. This line was the LP-0 status until 2026-09-07 and no check compared it with anything.
