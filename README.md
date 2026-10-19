# NonGravForces 0.16.0

Non-gravitational forces, torques and partials for LEO satellites (GRACE-FO reference case; classic GRACE-A/B macro model since N-07). Third concept, designed from scratch after the audits of `SatelliteNonGravToolbox` v1 and `NonGravToolbox` v0.1.1.

* Specification and plan: `HANDOVER.md` (single source of truth)
* Theory and implementation: `doc/src/NonGravForces_Theory_and_Implementation.tex` — the source is the deliverable. CI builds the PDF on every run and uploads it as the `theoryDocument` artefact; build it locally with `pdflatex` twice in `doc/src`. No PDF is committed: the one that used to be said "version 0.1.0" for three weeks because nothing compared it with anything (HANDOVER O-38). <!-- versionScan:allow -->
* Help browser: `setupNonGravForces`, then `doc NonGravForces`

```matlab
run('../referenceSystem/setupReferenceSystem.m'); run('../gracefoReading/setupGraceFOReading.m')
run('../odeSuite/gjToolboxSetup.m')      % optional for the toolbox, required by ciGate
setupNonGravForces; runNgfTests
```
Status: 59 catalogued units, all delivered, no stubs left (`auditPackage.py` compares this number with `doc/functionlist.csv`); all are verified in MATLAB. **Version 0.16.0** adds multi-day calibration with one common scale (N-09), measured temperature and magnetorquer regressors (N-11a), the calibration figures (N-05), a frozen regression baseline over January 2020 for GRACE-C and -D that every release must reproduce (N-10a), and a cheaper standard preset: the gas-stream shadow is read from a table (drag changes by at most 7e-13 m/s^2), which with a bitwise-neutral restructuring halves the model's share of a calibration. **437/437 green** on Tier A after the day bootstrap moved onto lsqvce's per-day contributions (R2026a Update 4, 2026-10-19; `tests/last_run.txt`).

Licensed MIT (see `LICENSE`). The NRLMSISE-00 C source is **not** in this repository — it carries no redistribution grant, so it is fetched from its maintainer instead:

```matlab
setupNonGravForces      % puts the toolbox on the path AND clones NRLMSISE-00 at the pinned commit
```

Setup is the only place in the toolbox that touches the network, it does so only when something is missing, and it prints a line when it does. It fetches the IERS and CelesTrak data files the same way. GRACE-FO Level-1B comes from the GFZ ISDC archive (147 MB for one day; setup keeps only the two files it needs). The toolbox is published via git only. `'FetchMsis', "never"` leaves the network alone; `NGF_NRLMSISE00` points at a copy you already have. `NonGravForces.buildMex` never fetches — it refuses and names setup. Everything except the `msis00` density backend works without any of it. See `THIRD-PARTY-NOTICES.md`. The authority is `doc/functionlist.csv`; `HANDOVER.md` Part A carries the level-by-level record.
