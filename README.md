# NonGravForces 0.9.1

Non-gravitational forces, torques and partials for LEO satellites (GRACE-FO reference case). Third concept, designed from scratch after the audits of `SatelliteNonGravToolbox` v1 and `NonGravToolbox` v0.1.1.

* Specification and plan: `HANDOVER.md` (single source of truth)
* Theory and implementation: `doc/NonGravForces_Theory_and_Implementation.pdf`
* Help browser: `setupNonGravForces`, then `doc NonGravForces`

```matlab
run('../referenceSystem/setupReferenceSystem.m'); run('../gracefoReading/setupGraceFOReading.m')
setupNonGravForces; runNgfTests
```
Status: 51 of 53 catalogued units delivered and verified in MATLAB, 2 documented stubs. **LP-7 is complete**; the four remaining stubs are blocked on data the project does not hold (a MAG1B file, a CERES granule, IGRF-14 coefficients). Next is LP-R, the release level.

Licensed MIT (see `LICENSE`). The NRLMSISE-00 C source is **not** in this repository — it carries no redistribution grant, so it is fetched from its maintainer instead:

```matlab
setupNonGravForces      % puts the toolbox on the path AND clones NRLMSISE-00 at the pinned commit
```

Setup is the only place in the toolbox that touches the network, it does so only when something is missing, and it prints a line when it does. It fetches the IERS and CelesTrak data files the same way. GRACE-FO Level-1B comes from the GFZ ISDC archive (147 MB for one day; setup keeps only the two files it needs). The toolbox is published via git only. `'FetchMsis', "never"` leaves the network alone; `NGF_NRLMSISE00` points at a copy you already have. `NonGravForces.buildMex` never fetches — it refuses and names setup. Everything except the `msis00` density backend works without any of it. See `THIRD-PARTY-NOTICES.md`. The authority is `doc/functionlist.csv`; `HANDOVER.md` Part A carries the level-by-level record.
