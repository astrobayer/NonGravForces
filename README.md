# NonGravForces 0.6.1

Non-gravitational forces, torques and partials for LEO satellites (GRACE-FO reference case). Third concept, designed from scratch after the audits of `SatelliteNonGravToolbox` v1 and `NonGravToolbox` v0.1.1.

* Specification and plan: `HANDOVER.md` (single source of truth)
* Theory and implementation: `doc/NonGravForces_Theory_and_Implementation.pdf`
* Help browser: `setupNonGravForces`, then `doc NonGravForces`

```matlab
run('../referenceSystem/setupReferenceSystem.m'); run('../gracefoReading/setupGraceFOReading.m')
setupNonGravForces; runNgfTests
```
Status: 43 of 53 catalogued units delivered and verified in MATLAB, 10 documented stubs; `kernel.singleEpoch` and the C kernel behind it (LP-6b) are delivered and measured against the frozen Python mirror, and await a MATLAB run report. The authority is `doc/functionlist.csv`; `HANDOVER.md` Part A carries the level-by-level record.
