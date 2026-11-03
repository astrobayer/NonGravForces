# HANDOVER — NonGravForces 0.18.0 (concept 3)

**Single source of truth** for the design, the audit that motivated it, the work plan, the rules, and the decisions of `NonGravForces`, a MATLAB toolbox computing non-gravitational forces, torques and partial derivatives for LEO satellites, with GRACE-FO accelerometer processing as the reference application and a sensor layer that admits future instruments (electrostatic, cold-atom, optomechanical, drag-free).

Revision 2.63 (release 0.18.0: the transient thermal walls in standard and precise, SOLAARS-CF corrected, the help audit's eight defects; PR #129 ☑ 511/511 on CI 41247178221) · 2026-10-31 · Claude (Anthropic) for Matthias Weigelt. Revision 2.62 (N-12, PR #124: IMU1B reader, ITSG attitude source, attitude-quality diagnostic; the gyro axes are not separable on a nominal arc) · 2026-10-30 · Claude (Anthropic) for Matthias Weigelt. Revision 2.61 (release 0.17.0: the observed density scale in the standard preset; baseline re-frozen and benchmark re-measured on the CI runner) · 2026-10-29 · Claude (Anthropic) for Matthias Weigelt. Revision 2.60 (PR #122 part 2: the standard preset takes drag.densityScale combined, strictly; release-checks workflow) · 2026-10-29 · Claude (Anthropic) for Matthias Weigelt. Revision 2.59 (PR #122: drag.densityScale on the native path -- a daily table in mex/ngfRhs.c; app.DynamicsSource no longer refuses it) · 2026-10-29 · Claude (Anthropic) for Matthias Weigelt. Revision 2.58 (PR #112: Swarm POD densities into DATAPOOL for an observed density scale, O-19; DTM2020 ☑ 459/459) · 2026-10-27 · Claude (Anthropic) for Matthias Weigelt. Revision 2.57 (PR #111: the DTM2020 month -- k_drag 0.96, 1.79 nm/s^2 along-track; the factor is the MSIS density, the gas-surface attribution corrected) · 2026-10-27 · Claude (Anthropic) for Matthias Weigelt. Revision 2.56 (PR #111: DTM2020 operational as density backend `dtm2020` -- CNES's unmodified subroutines fetched by setup, built by buildMex, refused on the native path, O-52) · 2026-10-27 · Claude (Anthropic) for Matthias Weigelt. Revision 2.55 (PR #111: the month -- ACT1B over all 31 days of 2020-01 C wants k_drag 0.44, GROOPS 0.48: the along-track drag is about twice too large) · 2026-10-27 · Claude (Anthropic) for Matthias Weigelt. Revision 2.54 (PR #109: drag diagnosis -- GROOPS's drag is 0.48 of this toolbox's, a constant scale; no single option reaches it) · 2026-10-27 · Claude (Anthropic) for Matthias Weigelt. Revision 2.53 (PR #108: io.spaceWeather reads Kyoto Dst and CLS F30, O-9 closed; O-19's comparison with TU Graz's GROOPS forces and a model-vs-ACT1B comparison after Huckfeldt et al. 2024, both on real days; gate 452) · 2026-10-26 · Claude (Anthropic) for Matthias Weigelt. Revision 2.52 (PR #107: observations of the reference days -- Kyoto final Dst, CLS F30, TU Graz attitude/non-conservative forces/orbits -- fetched by setupNonGravForces('FetchPool') into DATAPOOL and cut into committed excerpts on the CI runner; classic GRACE Level-1B parked, GFZ ISDC down, O-51) · 2026-10-26 · Claude (Anthropic) for Matthias Weigelt. Revision 2.51 (ACH1B prepared under its own name: timeSeriesPreparation >= 1.5; the label workaround removed, shown bitwise equivalent) · 2026-10-20 · Claude (Anthropic) for Matthias Weigelt. Revision 2.50 (N-10b: ACH1B as GRACE-D's oracle -- radial scale 0.784 -> 0.815, radial residual -15 %; N-10a's transplant attribution corrected; N-10 complete) · 2026-10-20 · Claude (Anthropic) for Matthias Weigelt. Revision 2.49 (N-11b: the temperature-dependent scale; regressor interpolation no longer bridges gaps; N-11 complete) · 2026-10-20 · Claude (Anthropic) for Matthias Weigelt. Revision 2.48 (fitCalibration's day bootstrap on lsqvce's per-day contributions, lsqvce >= 6.35: calibrate with 100 replicates 1410 -> ~340 s, identical sigmas) · 2026-10-19 · Claude (Anthropic) for Matthias Weigelt. Revision 2.47 (release 0.16.0: gate 435/435, benchmark re-measured, the regression baseline passed at 0.16.0; O-35 closed) · 2026-10-18 · Claude (Anthropic) for Matthias Weigelt. Revision 2.46 (the flow-shadow table's effect on a 10-day calibrate, measured: model half -46 %, bootstrap unchanged) · 2026-10-17 · Claude (Anthropic) for Matthias Weigelt. Revision 2.45 (geom.flowShadow Method "table", the standard preset's default, on both evaluation paths) · 2026-10-17 · Claude (Anthropic) for Matthias Weigelt. Revision 2.44 (occlusion and radiation force cheaper, every value bitwise unchanged: Model.accel -15 %) · 2026-10-16 · Claude (Anthropic) for Matthias Weigelt. Revision 2.43 (N-10a: where the baseline's four hours go, measured; the bootstrap claim corrected) · 2026-10-16 · Claude (Anthropic) for Matthias Weigelt. Revision 2.42 (N-10a: the regression baseline 2020-01 C/D frozen on the DATAPOOL; audit check 14) · 2026-10-16 · Claude (Anthropic) for Matthias Weigelt. Revision 2.41 (N-11a: measured regressors -- temperatures, magnetorquer currents -- in fitCalibration and Calibration; AHK1B/MAG1B fixtures) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.40 (N-09: multi-arc calibration -- one common scale, per-day local biases through lsqvce; three real days of GRACE-C fetched as fixtures) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.39 (N-05: `app.plotCalibration` and `app.plotAttribution`; attributeRevolutions' 20-revolution minimum documented, 'BlockLength' added) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.38 (release 0.15.0: O-49 closed; benchmark re-measured on Tier A) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.37 (N-03c: `app.perRevolution` and `app.attributeRevolutions`; N-03 complete; O-49 opened for release 0.15.0; the owner ratified every open decision of 2026-10-15) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.36 (N-03b: `app.compareCalibrations`, model comparison by blocked CV; `fitCalibration` gains 'Span' and info.prediction; an N-03a defect fixed) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.35 (N-03a: the time-varying scale in `app.fitCalibration` and `app.Calibration`; N-03 split into three PRs) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.34 (O-47 closed: the six deprecated `qc` forwarders deleted with their tests, catalogue rows and hushes; release 0.14.0, 54 units; benchmark reference re-measured) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.33 (release 0.13.0: N-02c verified 414/414 on the CI runner; `tools/benchmark.m` re-run on Tier A and its reference committed; O-48 closed) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.32 (N-02c and N-06: `app.calibrate` rewritten over `timeSeriesPreparation.run`; A.2 reworded; lsqvce and timeSeriesPreparation in the table; O-48 opened) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.31 (N-02b: `app.fitCalibration` on `lsqvce.fit`; lsqvce re-enters the dependency table; N-02a verified 408/408) · 2026-10-15 · Claude (Anthropic) for Matthias Weigelt. Revision 2.30 (N-02a: `app.syntheticPair` and `app.Calibration`; the audit's function-name parser fixed) · 2026-10-14 · Claude (Anthropic) for Matthias Weigelt. Revision 2.29 (N-01: `qc` forwards to timeSeriesPreparation; O-1 and O-18 moved, O-47 opened) · 2026-10-14 · Claude (Anthropic) for Matthias Weigelt. Revision 2.28 (N-04: Electrostatic's default kernel is the L1B CRN kernel from timeSeriesPreparation; timeSeriesPreparation and lsqvce are gate siblings) · 2026-10-14 · Claude (Anthropic) for Matthias Weigelt. Revision 2.27 (N-07: `geom.grace` delivered and verified on the CI runner, 388/388; O-46 opened) · 2026-09-22 · Claude (Anthropic) for Matthias Weigelt. Revision 2.26 (O-45 opened: the accelerometer calibration application stays in `app` and is upgraded in place, `qc` moves to the new shared toolbox timeSeriesPreparation, lsqvce re-enters; design in its concept) · 2026-09-18 · Claude Opus 5 (Anthropic). Revision 2.25 (io.igrfCoefficients delivered) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.24 (io.readMAG1B delivered; O-8 closed, O-37 opened) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.23 (mltbx deprecated; GRACE-FO L1B fetch corrected and wired; O-8 unblocked) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.22 (O-35 answered: optional data download in setup) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.21 (LP-R: WP-22 benchmark and packaging) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.20 (LP-R: WP-22 demos, run by the gate) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.19 (LP-R: WP-21 documentation build, O-11 closed) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.18 (LP-R opened: the C2 toolbox scan) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.17 (O-34: the fetch moves into setupNonGravForces, which owns the URL and the pin) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.16 (O-34 closed: NRLMSISE-00 fetched from upstream at a pinned commit, not redistributed) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.15 (LP-7 COMPLETE: sensor.DragFree delivered, MIT licence adopted, O-29 accepted; version 0.8.0) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.14 (LP-7: sensor.AtomInterferometer delivered; version 0.7.2) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.13 (owner ruling: lsqvce is overkill, WP-19 gets the toolbox's own solver; version 0.7.1) · 2026-09-09 · Claude Opus 5 (Anthropic). Revision 2.12 (LP-7: WP-29 sensor.Electrostatic delivered; version 0.7.0) · 2026-09-08 · Claude Opus 5 (Anthropic). Revision 2.11 (owner rulings recorded: LP-6b's two scope decisions approved, and P-11 — both evaluation paths are first class; version 0.6.3) · 2026-09-08 · Claude Opus 5 (Anthropic) for Matthias Weigelt. Revision 2.10 (LP-7 opened: `app.forceBudget` delivered, the gjToolbox and lsqvce surfaces read at last; version 0.6.2) · 2026-09-08 · Claude Opus 5 (Anthropic) for Matthias Weigelt · MATLAB reached through the CI runner, the MCP bridge still silent. Revision 2.9 (LP-6b single-epoch C/MEX delivered; version 0.6.1) · 2026-09-08 · Claude Opus 5 (Anthropic) for Matthias Weigelt · written on **Tier B (neither the MATLAB nor the filesystem MCP bridge answered)**: the C is measured against the frozen Python mirror by `mex/test/ngfKernelParity.c`, which needs no MATLAB, and the MATLAB side — `kernel.singleEpoch`, `tests/tNgfMex.m`, the second `buildMex` target — is ◐ until a run report is read. Revision 2.8 (LP-6 Model and partials delivered; version 0.6.0) · 2026-09-08 · Claude Opus 5 (Anthropic) for Matthias Weigelt · MATLAB R2026a Update 4. Revision 2.7 (O-12 discharged: the `l1bLoad` split removed, `dataEmpty` reachable for several files, O-25 opened; version 0.5.3) · 2026-09-08 · Claude Opus 5 (Anthropic) for Matthias Weigelt · written on **Tier B (no MATLAB reachable)**, then **verified on the self-hosted runner the same day: 245 passed, 0 failed, 0 incomplete, 0 audit shortfalls, predicted 245**. The ◐ it carried is discharged; 0.5.3 is ☑. Revision 2.6 (debt sweep: the category audit made honest, records reconciled with the tree; version 0.5.2) · 2026-09-07 · Claude Opus 5 (Anthropic) for Matthias Weigelt · MATLAB R2026a Update 4. Revision 2.5 (LP-5x, 0.5.1) · 2026-09-07. Revision 2.4 (LP-5, 0.5.0) · 2026-09-07. Revision 2.3 (LP-4, 0.4.0) · 2026-09-07. Revision 2.2 (LP-3, 0.3.0) · 2026-09-07. Revision 2.1 (LP-2, 0.2.0) · 2026-09-06. Revision 2.0 · 2026-08-27 · Claude Fable 5 (Anthropic).

Companion: `doc/src/NonGravForces_Theory_and_Implementation.tex` (theory, practical considerations, implementation), built by CI into the `theoryDocument` artefact and never committed as a PDF (O-38). This document does not repeat the physics; it cites the theory document by section.

---

## Part A — Read this first

### A.1 What was asked and what was delivered

Asked: scrutinise `SatelliteNonGravToolbox` (v1) and `NonGravToolbox` v0.1.1 in detail, design a third concept from scratch, sketch (not implement) the toolbox, include tests and documentation as part of the toolbox, produce a theory PDF, pay special attention to error handling, cover future sensors, and add outlier/impulse/jump detection.

Delivered (status key: ☑ verified in MATLAB, ◐ delivered-untested, ○ specified stub):

| Item | Where | Status |
|---|---|---|
| Audit of v1 and v0.1.1 (Part C) | this document | — |
| Architecture and workflow (Part D), function set (Part F) | this document, `doc/functionlist.csv`, `doc/html` | — |
| L0 core: `const`, `toolboxVersion`, `dependencyReport`, `fail`, `provenance`, `Options`, `Model` (constructor) | `+NonGravForces/` | ☑ |
| Geometry: `geom.MacroModel`, `geom.graceFO` (JPL D-56935 Tables 4–5) | `+NonGravForces/+geom/` | ☑ |
| Quality control: `qc.hampel` (with an independent standard-library oracle), `qc.eventTable` (schema) | `+NonGravForces/+qc/`, `tools/pyxcheck/` | ☑ |
| Sensor interface `sensor.Sensor` (abstract) | `+NonGravForces/+sensor/` | ☑ (interface) |
| Data I/O (LP-2, 2026-09-06): `io.readACC1B`, `readSCA1B`, `readGNV1B`, `readTHR1B`, `readMAS1B`, `massModel`, `spaceWeather`; real-data fixtures under `data/samples`, `data/spaceweather` | `+NonGravForces/+io/`, `data/`, `tests/tNgfIo.m` | ☑ (116/116, 2026-09-06) |
| Environment (LP-3, 2026-09-07): `env.density` (exponential, msis00 via the MEX of the vendored Brodowski C), `eclipseFactor` (cylindrical, conical, lunar), `earthRadiance` (pointSource, knocke cap tiles), `radiationSources`; `buildMex`; pure-Python oracles in `tools/pyxcheck` | `+NonGravForces/+env/`, `+NonGravForces/buildMex.m`, `mex/`, `tests/tNgfEnv.m` | ☑ (142/142, 2026-09-07) |
| Geometry and attitude (LP-4, 2026-09-07): `MacroModel` polygons and heat capacities, `graceFO` closed polygons and polygon centroids, `boxWing`, `selfShadow` (ray cast + icosphere LUT), `AttitudeProvider` (+ `nominal`), `checkAttitudeSense`, `thermalState`; `eclipseFactor` `solaarsCF` (O-13 closed) | `+NonGravForces/+geom/`, `tests/tNgfGeometry.m`, `tools/pyxcheck/selfshadow_check.py`, `solaars_check.py` | ☑ (169/169, 2026-09-07) |
| Kernel and quality control (LP-5, 2026-09-07): `gsiCoefficients` (sentman, maxwell, constantCd), `accommodation` (fixed), `dragPanelForce` (per species), `radiationPanelForce`, `thermalRecoil`, `assemble`, `bodyTorques`; `qc.residuals`, `detectImpulses`, `detectJumps` (PELT), `screen`; O-2 closed (unit normals) | `+NonGravForces/+kernel/`, `+NonGravForces/+qc/`, `tests/tNgfKernel.m`, `tests/tNgfQcScreen.m`, `tools/pyxcheck/sentman_check.py`, `pelt_check.py` | ☑ (231/231, 2026-09-07) |
| Model and partials (LP-6, 2026-09-08): the six `Model` evaluation methods, `env.wind` (co-rotation), `kernel.partials` (analytic per-panel derivatives) | `+NonGravForces/Model.m`, `+env/wind.m`, `+kernel/partials.m`, `tests/tNgfModel.m` | ☑ (262/262, 2026-09-08) |
| Debt sweep (2026-09-07) 244/244: the category audit reports shortfalls for exempted units too and three hidden coverage holes are filled; `qc.eventTable` refuses non-finite numbers; `auditPackage.py` checks the catalogue's field count and the vendored NRLMSISE-00 digests; O-10 closed upstream; the records reconciled with the tree | `tests/`, `auditPackage.py`, `HANDOVER.md`, `CLAUDE.md`, `README.md` | ☑ (244/244, 2026-09-07) |
| Cross-validation tooling (LP-5x, 2026-09-07): `tools/pyxcheck/ngfkernel/` independent standard-library mirror of the frozen kernel, parity driver and case list, `tests/tNgfMirror.m`; `hampel_check.py` made runnable | `tools/pyxcheck/`, `tests/tNgfMirror.m` | ☑ (238/238, 2026-09-07) |
| O-12 discharge (2026-09-08): `io/private/l1bLoad` routes every read through `gracefoReading.readL1bMulti`, one file or many; the `isscalar(files)` split and the two hand-fabricated plural aux fields are gone, `dataEmpty` is reachable for any number of files, the declared `gracefoReading` minimum moves to 1.0.1, O-25 opened | `+NonGravForces/+io/private/l1bLoad.m`, `dependencyReport.m`, `tests/tNgfIo.m` | ☑ (245/245, 2026-09-08; predicted 245) |
| Observations of the reference days (PR #107, 2026-10-26): `setupNonGravForces('FetchPool')` into DATAPOOL (Kyoto final Dst, CLS F30, TU Graz attitude / non-conservative forces / orbits), `tools/cutFixtures.py`, `.github/workflows/fetch-data.yml` on the self-hosted runner; excerpts in `data/spaceweather`, `data/samples/itsg`. No reader yet; classic L1B parked (O-51) | `setupNonGravForces.m`, `tools/`, `data/` | ☑ fetch and cut (run 40381190966); gate 444 unchanged |
| Readers and comparisons on the reference days (PR #108, 2026-10-26): `io.spaceWeather` families "dst" (Kyoto WDC) and "f30" (CLS, with the 81-day mean); `tNgfModel/modelGroopsReference` (Model against TU Graz's GROOPS-modelled forces on their own attitude and orbit, 2020-01-01 C: correlation 0.983 / 1.000 / 0.968, RMS ratio 1.50 / 1.18 / 1.20); `tNgfModel/modelAct1bReference` (Model against ACT1B, 3 x 30 min of 2020-01-01..03 C, one bias per day after Huckfeldt et al. 2024 eqs. 7-10: 5.5 / 3.4 / 2.1 nm/s^2 RMS) | `+NonGravForces/+io/spaceWeather.m`, `tests/tNgfIo.m`, `tests/tNgfModel.m` | ☑ (452/452, CI 40411254932) |
| DTM2020 density backend (PR #111, 2026-10-27): `env.density("dtm2020", ...)` and `Options` `drag.density = "dtm2020"` (opt-in, no preset), CNES's operational DTM2020 (F10.7, Kp) through `mex/ngfDtm2020Shim.F90` and `mex/ngfDtm2020Gateway.c`, built by `buildMex` target `dtm2020` from the files `setupNonGravForces` fetches (SWAMI MCM at a488a7c, SHA-256 pinned); `io.spaceWeather`'s celestrak query gains `f107Last81`, `kpDelayed3h`, `kpMean24h`; oracle CNES's four benchmark cases; `app.DynamicsSource` refuses it (O-52); `tools/dragScaleMonth.m` runs it beside MSIS | `+NonGravForces/+env/density.m`, `mex/ngfDtm2020*`, `tests/tNgfEnv.m` | ☑ (CI 40485043763: 459/459, predicted 459; the first run 40470246057 was 457/459 on two older contracts still asserting the 2026-09-15 refusal, fixed; month run 40470246020) |
| Swarm density scale (PR #112, 2026-10-27, O-19): ESA's Swarm A/B/C POD-only densities (`SW_OPER_DNSxPOD_2_`) fetched through VirES HAPI by `setupNonGravForces('FetchPool')` and cut into `data/samples/swarm` (first fetch-data run 40494585676: 2020-01, three satellites); new unit `io.readSwarmDensity` (header decoded, columns by name, units checked; 60 units); `tools/dragScaleMonth.m` gains three density-scale variants (Swarm B, A+C, interpolated in height to GRACE-FO), the preset's drag times rho_Swarm / rho_model per day | `+NonGravForces/+io/readSwarmDensity.m`, `tests/tNgfIo.m`, `tools/dragScaleMonth.m` | ☑ (CI 40505746004: 466/466, predicted 466; the first run 40500042251 was 466/466 passed but red on ONE audit shortfall -- the seven tests sat in the 'metamorphic' TestTags block, and the category is the block's tag, not the name; regrouped. Month run 40500042274: Swarm-interp k_drag 0.903, rms_x 2.08 nm/s^2) |
| Observed density scale as a Model option (PR #113, 2026-10-27, O-19): `Options` `drag.densityScale = "swarm"` (opt-in, no preset); new unit `env.densityScale` (per UTC day f = sum(rho_Swarm)/sum(rho_model) for Swarm A, B, C at their own positions; `.at(epoch, h)` interpolates log f linearly in height between A+C and B); `Model` 'DensityScale' multiplies density and every species density by `S.at(epoch, geodetic height)`, refuses a scale without the option or formed for another drag.density; `app.DynamicsSource` refuses the option (no native table); `tools/dragScaleMonth.m` runs the option as a Model variant | `+NonGravForces/+env/densityScale.m`, `Model.m`, `Options.m`, `tests/tNgfEnv.m`, `tests/tNgfModel.m` | ☑ CI 40528955139: 475/475, predicted 475 (written without MATLAB) |
| The density scale on the native path (PR #122, 2026-10-29, O-19, owner: standard takes "combined" and stays integrable): `app.DynamicsSource` no longer refuses `drag.densityScale`; `prepare` turns the Model's scale into one row per UTC day of the driver bins (`mex/ngfRhs.h` `scale`, 6 x days: kind 0 missing / 1 Swarm's bracket / 2 one factor, low, high, hLow, hHigh), and `mex/ngfRhs.c` multiplies rho and every species by it right after the density -- Swarm's bracket log-linear in height in env.densityScale's order, a driver or fallback day the factor `S.at` computed. The day is that of the ALIGNED driver bin, so the scale switches at the declared breakpoint of UTC midnight. A Swarm-only day without a factor is refused when reached (`rhsGateway:dataRange`), as Model.accel refuses it (`densityScale:dataDay`). Diagnostics carry `.densityScale` (NGF_RHS_DIAG 37). Three tests (parity vs Model.accel for swarm / driver / forced-fallback combined, the factor linearity, the refusal) | `+NonGravForces/+app/DynamicsSource.m`, `mex/ngfRhs.c`, `mex/ngfRhs.h`, `mex/ngfRhsGateway.c`, `tests/tNgfIntegration.m` | ☑ CI 40880333732: 493/493, predicted 493 (written without MATLAB; C syntax-checked with gcc) |
| The standard preset takes the combined density scale (PR #122 part 2, 2026-10-29, O-19, owner ruling): `Options` standard `drag.densityScale = "combined"` (fast, precise "none"); strict -- a standard Model evaluates only with 'DensityScale'; `NgfTestCase.standardUnscaled()` for tests about something else, `tNgfIntegration` on the preset as shipped; demos, benchmark and regressionBaseline build the scale (excerpts / pool); dragScaleMonth's reference row stays unscaled; `release-checks.yml` runs baseline and benchmark on the runner | `+NonGravForces/Options.m`, `tests/*.m`, `demos/demoForceBudget.m`, `tools/benchmark.m`, `tools/regressionBaseline.m`, `tools/dragScaleMonth.m`, `.github/workflows/release-checks.yml` | ☑ CI 40896217330: 493/493, predicted 493 (first run 40884016269 492/493: the gravity-free arc of DynamicsSourceMetamorphic left the Swarm bracket; it runs on the driver scale now) |
| N-12 attitude sources and the attitude-quality diagnostic (PR #124, 2026-10-30, #69): IMU1B fetched into the pool months and cut (30 min of 2020-01-01..03 C); `io.readIMU1B` (FiltAng = each gyro's filtered accumulated angle, deg -> rad, 8 Hz per gyro, gyro ids kept apart, GPS time enforced); `geom.AttitudeProvider.fromItsg` (TU Graz's attitude conjugated once); `geom.attitudeQuality` (SCA1B rates vs gyro angle differences on the same 1 s intervals, fitted rate = a.omega + b, residual by exact projection, nominal-frame angles with 'Orbit', 1/2 cpr); oracle `tools/pyxcheck/attitude_quality_check.py`. Finding: on an Earth-pointing arc the gyro axis, scale and bias are not separable; the residual (3.6-4.3e-7 rad/s RMS) is the result | `+NonGravForces/+io/readIMU1B.m`, `+geom/attitudeQuality.m`, `+geom/AttitudeProvider.m`, `tests/tNgfAttitude.m`, `tools/cutFixtures.py`, `setupNonGravForces.m`, `tools/pyxcheck/attitude_quality_check.py` | ☑ CI 40980858655: 509/509, predicted 509 (493 + 16); first run 40978582699 507/509 (oracle quaternions unnormalised; an epoch asserted below the product's 1 us) |
| Transient thermal on the native path (PR #126, 2026-10-31, owner's ruling: standard goes transient and stays integrable): `app.DynamicsSource` 'ThermalOrbit' -- prepare integrates the wall temperatures once along the reference orbit with Model.accel and pchips them onto the table nodes; `mex/ngfRhs.c` NGF_THERMAL_TABLE reads them by time and forms the recoil; refused without the orbit, with 'Native' false, or when the orbit does not cover the tables | `+app/DynamicsSource.m`, `mex/ngfRhs.c`, `mex/ngfRhs.h`, `mex/ngfRhsGateway.c`, `tests/tNgfIntegration.m` | ☑ CI 41154819317: 510/510, predicted 510 |
| The standard preset goes transient (PR #126, 2026-10-31, owner's ruling): `Options` standard thermal.model "transient" with Woeske's capacities (fitted, not improved: run 41140282748); `Model.accel` 'T0_K'; `forceBudget`/`densityInversion` carry the walls across chunks; tests: `standardUnscaled` steady, `tNgfIntegration` "standardSteady" for 'Native' false, standard as shipped through 'ThermalOrbit' | `+NonGravForces/Options.m`, `Model.m`, `+app/forceBudget.m`, `+app/densityInversion.m`, `+geom/graceFO.m`, `tools/benchmark.m`, `tests/NgfTestCase.m`, `tests/tNgfIntegration.m`, `tests/tNgfModel.m` | ☑ CI 41154819317: 510/510, predicted 510; first run 41152904075 509/510 (the between-node bound of the time tables met the temperature table, recorded apart) |
| User guide follows the transient standard preset (PR #127, 2026-10-31): preset table, the two owner-ruled exceptions to opt-in, the native-path refusal now only without 'ThermalOrbit' or with 'Native' false, a TRP paragraph on the measured effect, the Cost table marked as predating the transient walls (re-measured at the next release); the PDF is not rebuilt (O-38) | `doc/src/NonGravForces_UserGuide.tex` | documentation only, gate unchanged at 510 |
| Documentation pass (PR #128, owner 2026-10-31): help-vs-code audit of all 64 units (comment lines only); theory document updated (density scale section, DTM2020, transient walls, IMU attitude check, workflow, rsApiRotation claim removed); generated test catalogue `doc/src/tests_generated.tex` as the theory PDF's appendix, compared by audit check 6; GettingStarted.html rewritten; `tools/userGuideTranscript.m` + release-checks `transcript`. CODE FIX found by the audit: transient walls held the interval across every chunk boundary in `forceBudget`, `densityInversion` and `Model.series`; later chunks now start one epoch early and drop it | `+NonGravForces/Model.m`, `+app/forceBudget.m`, `+app/densityInversion.m`, 60 help blocks, `tools/buildDocs.py`, `auditPackage.py`, `doc/src/*`, `tests/tNgfApp.m` | ☑ CI on 8734694: 511/511, predicted 511 (510 + `tNgfApp/forceBudgetTransientChunks`), 0 audit shortfalls |
| The eight code defects of the PR #128 help audit (PR #129, owner 2026-10-31): `env.densityScale` refuses TOLEOS files (were silently dropped), `env.densityDriverScale` refuses the unadjusted F30, `io.igrfCoefficients` refuses a multi-epoch struct, `provenance` hashes with exact FNV-1a (the double product rounded; hashes change), `geom.selfShadow` reports a reused table's samples, `geom.flowShadow` table fills `info` on every path, `AttitudeProvider.fromItsg` sorts and de-duplicates (refuses conflicting copies), `io.readMAG1B` records its files | `+env/densityScale.m`, `+env/densityDriverScale.m`, `+io/igrfCoefficients.m`, `provenance.m`, `+geom/selfShadow.m`, `+geom/flowShadow.m`, `+geom/AttitudeProvider.m`, `+io/readMAG1B.m`, tests in tNgfEnv/Io/Core/Geometry/Attitude | ☑ CI 41245260551: 511/511, predicted 511 (assertions inside existing methods) |
| The precise preset goes transient (PR #129, 2026-10-31, owner's ruling on the recommendation of run 41213060789): `Options` precise thermal.model "transient" (steady state led ACT1B by 12-20 s on y/z, transient 2-7 s); no baseline is frozen for precise; the native path serves it through 'ThermalOrbit', which `tNgfIntegration` already passes for every transient configuration (precise, walker2014, sesam now among them). **standard keeps srp.eclipse "conical"** (owner: leave as is, no re-freeze). Help (Options, thermalRecoil, thermalState), theory and user guide follow; `tNgfCore/optionsContract` asserts the thermal ladder | `Options.m`, docs, `tests/tNgfCore.m` | ☑ CI 41247178221: 511/511, predicted 511, Code Analyzer clean |
| 10 documented stubs raising `NonGravForces:<name>:notImplemented` (**superseded: none left since LP-C; this table is the level record up to LP-6b, the status lines in README and CLAUDE.md section 5 are current**) | all packages | ○ |
| Test infrastructure: `NgfTestCase`, `runNgfTests`, `auditTestCategories`, exemptions; 11 test classes, 262 tests, **gate GREEN 262/262** (74 at LP-0, 116 at LP-2, 142 at LP-3, 169 at LP-4, 231 at LP-5, 238 at LP-5x, 244 after the debt sweep) | `tests/` | ☑ |
| CI: `ci.yml`, `ciGate.m`, `auditPackage.py`; Code Analyzer clean (0 findings) | repository root, `.github/workflows/` | ☑ (local run 2026-08-27) |
| Help-browser documentation (`info.xml`, `helptoc.xml`, six pages) | `doc/html/` | ☑ (`doc NonGravForces` opens) |
| Theory & implementation | `doc/src/*.tex` (the deliverable; CI builds the PDF as an artefact, `auditPackage.py` refuses a committed one) | — |

Run log of the verifying run: `tests/last_run.txt` (2026-09-08, R2026a Update 4, 262 passed, 0 failed, 0 incomplete, 0 audit shortfalls, predicted 263 then 262 after merging main — see Part I; debt sweep: 244; LP-5x: 238; LP-5: 231; LP-4: 169; LP-3: 142; LP-2: 116; LP-0: 74).

### A.2 The rule that governs every other rule

> **NonGravForces computes forces and torques, and applies them; it implements nothing it can call.**

*Reworded 2026-10-15 under O-45 N-06 (#63). Until then the rule read "computes forces and torques and nothing else", and the `app` layer had outgrown it: `app.calibrate` calibrates the accelerometer against those forces (O-45). The application stays here because the a priori is this toolbox's; what it may not do is re-implement a neighbour. Time and frames are referenceSystem's, reading is gracefoReading's, the instrument stages are timeSeriesPreparation's, integration is gjToolbox's, and every least-squares number is lsqvce's.*

Five toolboxes of the same owner form the ecosystem; each is called by qualified name and none is copied, vendored, wrapped or aliased (versions as in `NonGravForces.dependencyReport`, the authority):

| Toolbox | Owns | NonGravForces uses |
|---|---|---|
| `referenceSystem` ≥ 2.3 (raised from 2.1 by N-02c) | time scales, EOP, ITRS↔GCRS with derivatives, geodetic coordinates, quaternion algebra, attitude interpolation/integration, Sun/Moon/planet states (SOFA analytic or DE440 via MICE), int64 GPS nanoseconds | everything kinematic; `app.calibrate` converts Arc time with `epochTo(..., "gpsNanoseconds")` |
| `gracefoReading` ≥ 1.1 (1.0.1 until 2026-10-26, raised for `readGroopsInstrument` on real GROOPS output) | L1B ASCII reading, quality flags, GROOPS instrument export/import | all L1B input; TU Graz's GROOPS products in `tNgfModel/modelGroopsReference` (O-19) |
| `timeSeriesPreparation` ≥ 1.5 (optional for the toolbox, required by `ciGate`; 1.4 until 2026-10-20, raised for the ACH1B preset) | the `Arc`, the CRN kernel and its application, mask, align, QC (`qc` moved there, N-01), noise, observability, bootstrap, `run` | `sensor.Electrostatic`'s default filter (N-04); `app.calibrate` runs the instrument stages through `run`; `app.fitCalibration` filters every column with `crnApply` and gates with `observability` and `bootstrap` |
| `lsqvce` ≥ 6.35 (optional for the toolbox, required by `ciGate`; 6.22 until 2026-10-19, raised for `R.contributions`) | batch least squares (declaration kinds A/B/C), VCE, ARMA/PSD whitening, robust reweighting, diagnostics | `app.fitCalibration` declares kind A per axis; its day bootstrap solves the fit's per-day contributions (`lsqvce.solveContributions`), its single-arc block bootstrap refits with kind B; `app.calibrate`'s attribution solves with it. **The row is true again since N-02b**: the toolbox's own QR of 0.12.0 is retired (B.4) |
| `gjToolbox` (odeSuite) ≥ 0.1 (optional for the toolbox, required by `ciGate`) | Gauss-Jackson orbit integration; `gjToolbox.DynamicsSource` interface | `app.DynamicsSource` serves it the non-gravitational acceleration (`Model.rhs` is gone, O-4) |

Consequences that differ from v0.1.1: **no gravity, no equations of motion, no integrator, no least-squares code** inside NonGravForces (v0.1.1 planned WP-17 `dyn/gravity` + `eomFcn` and hand-rolled calibration; both are removed, see C.3). The `app` layer composes; it does not solve.

### A.3 How to run

```matlab
run('C:\Users\matth\Documents\MATLAB\referenceSystem\setupReferenceSystem.m')
run('C:\Users\matth\Documents\MATLAB\gracefoReading\setupGraceFOReading.m')
cd('C:\Users\matth\Documents\MATLAB\NonGravForces'); setupNonGravForces   % adds root, tests/, tools/; registers doc
runNgfTests                          % gate; 'ExcludeSlow',false adds the speed tag
doc NonGravForces                    % help browser
```
`setupNonGravForces` adds `tests/` itself (v0.1.1 did not: finding V-13 friction, closed by design).

**Continuous integration** (`.github/workflows/ci.yml`, modelled on referenceSystem's): job `audit` runs `auditPackage.py` (stdlib only: catalogue ↔ files, help template, `fail` taxonomy, banned constructs); job `matlab` checks out `referenceSystem` and `gracefoReading` as siblings, runs `codeIssues` on the package (must be empty) and then `ciGate`, the one entry point that installs, runs the suite and judges. `ciGate` is also the local postflight (`matlab -batch ciGate`). Branch `main` is protected by a **ruleset**: both checks required (strict), `deletion` and `non_fast_forward`, no bypass actors. Classic branch protection is NOT in use, so the `/branches/main/protection` API answers "Branch not protected" — an instrument that cannot see rulesets, and reading only it would wrongly report `main` as unprotected. **No rule requires a pull request**; changes arriving by PR is a convention here, not a guardrail. A second ruleset protects `~ALL` tags against deletion and non-fast-forward.

**When the MCP bridge is down, CI is still a MATLAB tier** (learned 2026-09-08, LP-6b): the self-hosted runner answered a whole session in which neither the MATLAB nor the filesystem bridge did. Push the branch, let the `matlab` job run `ciGate`, and read `tests/last_run.txt` from the run artefact. It is slower than the bridge and it cannot be used to explore, but it turns a Tier B session's ◐ into a read run report, which is what P-1 asks for.

### A.4 Session protocol (unchanged in spirit from v0.1.1, tightened)

1. One level = one session = one deliverable zip. Paste the standing instruction (E.1), the level prompt (Part G) in full, and this document.
2. Nothing is ☑ without a MATLAB run report read by the author. `runNgfTests` writes `tests/last_run.txt`; the tick rests on that file, not on the pass count.
3. Predict the test count before the run; reconcile after (P-2).
4. A test failure is decomposed before any bound moves (R6). In this very session two failures were test-authoring errors (a typed oracle sum, a 0-based index) and one was a bound set by guess (interpreted vs compiled speed ratio); each is recorded in the test's comment with the measured number.

---

## Part B — Dependency contracts (measured 2026-08-27)

### B.1 referenceSystem 2.1.216 — functions NonGravForces calls

| Need | Call | Convention |
|---|---|---|
| epoch from L1B | `daySecondToEpoch(Nx2,"GPS")` | epoch struct `.day .sec .scale .source` |
| scale conversion | `convertTimeScale(ep,"TT"|"UT1"|"UTC")` | UT1 needs EOP |
| ITRS→GCRS | `rotationItrsToGcrs(ep,'Eop',eop)` → `[R,Rdot,Rddot]` | 3×3×N pages, ITRS→GCRS, tidal corrections default on, samples EOP beyond the arc (±3 d margin). **Rows before 1972-01-01 are trimmed with a warning** (re-measured 2026-09-07 against 2.1.216, O-10 closed): the shipped `eopc04_20.1962-now` now reads whole, 3652 of 23562 rows ignored, and `readGNV1B` works with it. Until 2.1.216 it raised `referenceSystem:leapSecondTable:beforeUtcEpoch` instead. `finals2000A.all` and the trimmed C04 agree to 9 mm / 5e-6 m/s on the 2020-01-01 arc. Cost measured: 600 epochs with `Rdot` 0.14 s, without 0.024 s |
| apply | `applyRotation(R,Rdot,r,v)` | carries ω×r; never formed by hand |
| geodetic | `cartesianToGeodetic(Nx3,'Ellipsoid',"WGS84")` | **rows in, default GRS80** → WGS84 passed explicitly |
| Sun/Moon | `geocentricBodyState("sun",ep)` / `Source="mice"` | ICRS-oriented 3×N; analytic tier RMS 8.3e3 m; DE440 tier via furnished kernels (de440s.bsp present) |
| attitude | `interpolateQuaternions(4×M,nodeEp,queryEp,'Method',"slerp"|"squad")`, `rotationFromQuaternion`, `multiplyQuaternions`, `invertQuaternion`, `attitudeConvention()` | scalar-first Hamilton, passive; SCA1B is ICRF→SRF, NonGravForces stores BODY→GCRS, conjugated once in `geom.AttitudeProvider` |
| angular rate | `angularVelocityFromRotationRate(R,Rdot)` | CoM-coupling term in `sensor.Electrostatic` |
| orbit frame | `rotationGcrsToRsw(r,v)` | budgets and plots |
| constants | `iersConstants()` (`.speedOfLight .earthRotationRate .gmEarth …`), `ellipsoidParameters` | no `au`: NonGravForces owns it in `const` |
| version | `toolboxVersion().version` | **struct, not string** — the first probe of this session was wrong and was corrected by measurement |
| C API | `mex/rsApi.h` `rsApiRotation` (0.56 µs) | **NOT linked, and the plan that said it would be was wrong** (read 2026-09-08): `rsApiRotation` is ITRS→GCRS, and the frozen kernel is body frame throughout, so the single-epoch MEX has no rotation to evaluate. The header stays the right route for a C *integrator* that needs R(t) in its own right-hand side; that is a caller of this toolbox, not this toolbox |

### B.2 gracefoReading 1.0.0 (measured 2026-09-06; was 2.3.0 on 2026-08-27, renumbered by its owner)

`readL1bMulti(files)` → `[time, obs, aux]` with `aux.products`, `aux.fileNames`, `aux.recordsPerFile` (no `product` field); `checkQuaternionNorm(obs, aux)` → `[ok Nx1, stats]`, default tolerance 1e-6, measured 9e-11 on the real SCA1B sample; `readL1b(file)` → `[time Nx2 GPS, obs, aux]`; `readL1bMulti`; `decodeQualFlags`; `checkQuaternionNorm`; `findGaps`; `version()`. Rules learned in v0.1.1 and kept: `obs` holds measurements only (flags, char columns, `_id` columns, time tags go to `aux`); `prod_flag` products have per-record columns (bit 0 rightmost); **MAS1B is cold-gas tank mass (~30 kg), not spacecraft mass (601.214 / 601.209 kg)** — misuse gives 20–40× accelerations; `io.massModel` defaults to the handbook constant and asserts > 500 kg.

### B.3 gjToolbox (odeSuite) — **superseded 2026-09-16: the interface exists and is served**

**Current state (read at `6a3f382`, used by LP-C WP-C9).** odeSuite delivers `gjToolbox.DynamicsSource` (construct / `prepare(context)` / `accel(t, R, V)`, `StateSize`, `ParameterSpec`, `Breakpoints`, `NativeEntry`), `gjToolbox.Dynamics` with `Binding` "vectorized" or "native", `gjToolbox.propagate`, `gjToolbox.Grid`, and the C contract `csrc/gj_native.h` (ABI 1: `gj_native_entry`, `gj_native_call`, `GJ_NATIVE_OK`/`EDOMAIN`/`EFAIL`). `NonGravForces.app.DynamicsSource` subclasses it and publishes a native entry; `Model.rhs` is deleted. gjToolbox is optional for this toolbox and **required by `ciGate`**, which runs `odeSuite/gjToolboxSetup.m` from the sibling checkout — a gate that skipped the integrator interface would be green about nothing.

**The 2026-09-08 record, kept because its lesson outlived its facts.** `astrobayer/odeSuite` was then at **level L0**: version metadata, validators, `Budget`, `Phase`, `Rat`, the coefficient tables and the build system. **There is no integrator and no dynamics interface.** Its README says the dynamics interface is level L3 and calls it `gjToolbox.DynamicsSource`; its own validator `mustBeDynamics` calls it `gjToolbox.Dynamics`, and the two disagree inside that toolbox. The live check is `mustBeDynamics`, which **degrades on purpose**: with no such class on the path it accepts any scalar object or function handle, and its documented example is `@(t,r,v)` — **three arguments, numeric t**.

`Model.rhs()` returns `@(epoch, r, v, q)`: **four** arguments, the first a `referenceSystem` epoch struct and the last a quaternion. It would pass `mustBeDynamics` today (it is a function handle) and then be called with the wrong signature by whatever L3 becomes. The record here said "`Model.rhs()` returns `@(t,r,v) a_ng`" until 2026-09-08; the code and this document had disagreed since LP-6 and nothing compared them (O-4, now narrowed with evidence).

What that meant for LP-7 and beyond: an adapter was needed, not a change to `Model` — something closing over the epoch origin, the attitude provider and the EOP. **That is exactly what shipped**, as a class rather than a handle, and `Model` was not touched. The lesson that survives: this record claimed for four months that `Model.rhs()` returned `@(t,r,v)`, the code had said otherwise since LP-6, and nothing compared the two. A dependency's surface — and one's own — is a guess until it is read.

### B.4 lsqvce — a dependency again since N-02b (2026-10-15); the 2026-09-09 removal below is history

**Re-entry, recorded under N-06 (#63).** The owner's concept for accelerometer calibration (O-45; `timeSeriesPreparation` `docs/CONCEPT.md` Q2) makes `lsqvce.fit` the one estimator every application calls, and N-02b met the re-entry condition the removal wrote: variance components, spectral noise models, robust reweighting and streaming accumulation now have work to do (the CRN-filtered residual is coloured; multi-arc calibration, N-09, accumulates arcs). `app.fitCalibration` declares kind A (one system, one day, one arc per contiguous run) per axis; its bootstrap refits kind B. The surface was read from source at 6.22.0 (`fit`, `Options`, `Results`, the API reference chapters `01_data` to `06_results`), which is the minimum. `app.calibrate`'s own QR solver of 0.12.0 is retired with the matrix interface (N-02c).

**The removal of 2026-09-09, kept as written:**

**`lsqvce` is overkill for this estimation problem and is gone from the dependency table.** The record here spent two revisions being wrong about that toolbox's surface — it named the 1.0 flat API long after 5.3 had replaced it with five user names — and the correct response was not a better record but the question of whether the dependency was needed at all. It is not.

Why not, stated so it is not re-litigated: `app.calibrate` fits `y = bias(t) + S a_model`, which is LINEAR in every unknown and **separates by axis** — row *i* of the measurement depends on row *i* of the bias coefficients and row *i* of `S` and on nothing else. That is three independent least squares of a few dozen unknowns, solved by QR on an equilibrated design matrix. Variance components, spectral noise models, streaming accumulation, arc factors and Helmert scaling — the things lsqvce exists for — have no work to do here.

What was given up, recorded so nobody discovers it as a surprise: the formal sigmas are computed under a diagonal weight matrix, and accelerometer noise is not white. `app.calibrate` says so in its own help and points at `qc.residuals` for whitening. If a future unit genuinely needs variance components, lsqvce comes back into the table; until then a table naming a toolbox nothing calls would be believed, and that is worse than an absent row.

### C.1 SatelliteNonGravToolbox v1 — verified defects (numbers re-derived, see PDF §2)

| # | Defect | Evidence in source | Magnitude |
|---|---|---|---|
| D1 | Sentman: exponential group of C_p and the whole shear term divided by s instead of s² | `satDragAccel.m` lines `Cp = ... / s_ratio`, `Ctau = ... / (s_ratio*sqrt(pi))` | shear ×8 at s≈8 (LEO). **Annotated 2026-09-07 (LP-5x):** this row describes v1's code in the drag/lift convention it used. In THIS toolbox's normal/tangential convention the tangential coefficient is l[E/(s√π) + γZ] and has no 1/s² term at all; the Z/(2s²) belongs to the lift coefficient. Three independent derivations agree, and the edge-on limit settles it |
| D1b | ambient gas treated as air, R = 287 J/kg/K | `R_specific = 287` | speed ratio wrong by √(29/16) = 1.35 at 500 km (atomic O); Sentman coefficients wrong by 10–20 % |
| D2 | spherical-Earth geodetic conversion | `eciToGeodetic.m` | polar altitude −21.4 km ⇒ density ×1.40 |
| D3 | two solar constants (SRP implied 1367, Earth-radiation 1361) | `satSRPAccel.m`, `satEarthRadiationAccel.m` | 0.44 % inconsistency |
| D4 | GPS tags used as UT1; mean-of-date Sun mixed with GCRF; no ITRS→GCRS path (GNV1B is Earth-fixed) | `jd_ut1` arguments everywhere, `sunPositionECI`, no EOP | 18 s Earth rotation (0.075°), 0.36° precession |
| D5 | NRLMSISE-00 fed adjusted and same-day F10.7 | `getSpaceWeather.m` | ±3.4 % seasonal |
| D6 | silent fallback exponential density when NRLMSISE absent (warning only) | `atmosDensity.m` | 10× fidelity change unannounced |
| D7 | the drag routine computes the Sun position itself to articulate arrays | `satDragAccel.m` → `sunPositionECI` | hidden coupling; two Sun computations per epoch |
| D8 | optical flat-plate model in five divergent copies; torque path lost the vis/IR split | `panelForcesECI`, `satSRPAccel`, `satEarthRadiationAccel`, `computeNonGravTorque`, `nongrav_all_mex.c` | 8 % torque error, every force evaluated twice |
| D9 | single-epoch (3×1) interfaces only; no vectorisation over epochs | all `src/core` | series cost = K × MATLAB call overhead |
| D10 | hand-rolled test harness (`assertClose`, PASS/FAIL prints), no `matlab.unittest`, no tags, no categories, "I could not execute the code" in the install guide | `tests/` | untestable claims |
| D11 | duplicated quaternion/time helpers (`quat2dcm_local`, `gmstFromJD`) instead of a reference-system dependency | `src/utils` | second authority for frames |

What v1 got right and is kept: the GRACE-FO macro model with quoted JPL D-56935 values and energy-closure check; the vis/IR split; tiled Earth radiation with self-shadowing clusters; dual-cone penumbra; MAS1B scholarship; physics-based attitude-sense check; caveat culture.

### C.2 NonGravToolbox v0.1.1 — design flaws found in this audit

v0.1.1 fixed D1–D11 in its specification and is well engineered where it is implemented (151 green tests, strict conventions). The flaws are architectural and organisational:

| # | Flaw | Where | Why it matters | Decision |
|---|---|---|---|---|
| F1 | **Scope creep: gravity and equations of motion inside a non-gravitational toolbox** (WP-17 `dyn/gravity`, `eomFcn`; WP-17b full RHS in C) | HANDOVER Part F LP-6/6b, `+dyn` | a second authority for gravity next to the owner's own shAnalysis/gjToolbox; the integrator toolbox forbids exactly this (its R9) | removed; `Model.rhs()` is the seam (A.2) |
| F2 | **Compatibility layer to v1 (WP-20 `compat/`, parity baselines)** | LP-R | parity with a toolbox whose numbers are wrong by ×8 is not a feature; it costs a session and freezes defective signatures | removed; v1 defects become regression tests with the *correct* numbers |
| F3 | **`env.sunGCRF` is a wrapper in disguise** (a "selector" over `referenceSystem.geocentricBodyState`) | `+env/sunGCRF.m` (311 lines) | violates A.2 of its own handover: one more authority for the Sun | removed; the tier is `Options.sun.source`, the single call site is `env.radiationSources` |
| F4 | **The handover mixes specification with run logs** (4171 lines; Parts L and K are session diaries) | HANDOVER.md | a reader cannot find the spec; the diary grows without bound | this document is specification + decisions only; run logs are `tests/last_run.txt` and git history |
| F5 | `setupNGTToolbox` does not add `tests/`; runner is 1327 lines | `setupNGTToolbox.m`, `runNgtTests.m` | recurring V-13 friction; the runner carries reporting that belongs to `matlab.unittest` plugins | setup adds `tests/`; runner ≤ 100 lines, gate = failures + audit + written log |
| F6 | Only four test categories; no robustness, metamorphic or speed category as a *requirement* | E.4 | NaN/gap behaviour and invariants were tested ad hoc | seven categories, counted mechanically |
| F7 | TSI oracle 1367 W/m² (Behzadpour Table 1) while `const.S0` = 1361 | memory/decision log vs `const.m` | the same D3 class of inconsistency, one level up | S0 = 1361 (IAU 2015 B3, SORCE/TSIS); 1367 is an `Options.srp.tsi` override for legacy comparisons, never a constant |
| F8 | `fast` preset uses a "Vallado + precession" Sun although the dependency gives the analytic SOFA Sun for free | D.4 | a second Sun implementation | removed |
| F9 | Empty placeholder packages `+tf`, `+dyn`, `+app`, `+geo`, `+kernel` with `Contents.m` only | tree | clutter without contracts | every package ships its catalogue row and stubs that raise |
| F10 | No sensor layer: calibration written directly against ACC1B | WP-19 | a second instrument (CAI, MAGIC) would duplicate the calibration | `sensor.Sensor` interface; calibration and QC coded against it |
| F11 | No quality control: outliers/impulses/jumps left to the calibration residuals | — | jumps absorbed by bias splines silently; impulses inflate σ₀ and mask blunders (lsqvce `snoop` documents exactly this masking) | `qc` package (LP-5) |
| F12 | Error handling by convention only (identifier format), no mnemonic taxonomy, no correction hint, no single raise point | E.3 style rule | tests cannot distinguish contract from data from numeric failures; users get messages without a fix | `NonGravForces.fail` with taxonomy and mandatory ` -> ` hint (Part E) |
| F13 | Density models planned as one class per backend (`DensityMSIS`, `DensityDTM`, `DensityMSIS2`) with a base-MATLAB NRLMSISE-00 port | WP-05/25 | three classes for one contract; a MATLAB port of MSIS is the slowest possible inner loop | one dispatcher `env.density(model, …)` with backends as private functions; MEX for the density inner loop is the sanctioned C use (C3) |
| F14 | Self-shadow LUT over an icosphere with 642 nodes fixed | WP-10 | resolution not tied to an accuracy claim | `Options.geom.lutResolution` with a convergence test in LP-4 |

Kept from v0.1.1 without change: the dependency contract discipline (Part B), the P1–P7 principles, the five failure modes, the `arguments`-block rule, the help template, the ◐/☑ discipline, the Python cross-check idea (`tools/pyxcheck`), the `obs/aux` rule, the MAS1B rule, WP numbering where the content survived.

### C.3 Approximations and deviations from best methods, and whether they are justified

| Approximation | Justified? | Toolbox position |
|---|---|---|
| 8-plate macro model instead of finite-element ray tracing | yes for calibration at the 10 % force level (Wöske 2019); not for gravity-quality SRP | default; `geom.selfShadow` raycast + future kind "mesh" hook |
| Sentman/DRIA with constant α | yes at 400–500 km (March 2019 α≈0.85); SESAM for solar-cycle dependence | `drag.accommodation = fixed | sesam | langmuir` |
| co-rotating atmosphere, no wind | cross-track only; HWM14 opt-in | `drag.wind` |
| TRP folded into SRP (2/3 factor) | crude; TRP is the second-largest radiation term (Siemes 2023) | `thermal.model = foldedIntoSrp | steadyState | transient` |
| conical penumbra, spherical Earth | good enough for drag work; SOLAARS-CF for SRP at the nm/s² level and for integrator smoothness | `srp.eclipse` |
| Knocke zonal albedo/IR | fine to ~0.2 nm/s²; CERES for anything better | `erp.backend` |
| constant handbook mass | standard practice; MAS1B-driven variable mass opt-in | `io.massModel` |
| CoM-proof-mass offset ignored | ≤100 µm × ω̇ ~ 1e-6 rad/s² ≈ 1e-10 m/s² — at the noise floor; included in the electrostatic sensor model, not in the physics | `sensor.Electrostatic` |
| antenna thrust, Lorentz force | ~1e-11 m/s² or below; hooks only, off by default | `Options` (LP-3b) |

### C.4 Limitations of the method (what no toolbox fixes)

Thermospheric density uncertainty of 10–30 % (model spread; NRLMSIS 2.1 ~10 % low vs MSISE-00, Guo & Zhang 2024) dominates along-track; cross-track is dominated by winds (HWM14 RMS 40–80 m/s) and TRP; the penumbra transition and eclipse thermal transients cannot be modelled below ~1 nm/s² without a transient thermal model; accelerometer artefacts (twangs, heater switching, magnetorquer spikes) are not physics and are handled by `qc`, not by the model.

---

## Part D — Architecture

### D.1 Principles (P1–P7 kept from v0.1.1, P8–P10 new)

P1 one authority per fact · P2 one kernel (per-panel forces; force and torque are reductions) · P3 sources, not effects · P4 body frame inside · P5 configure once, evaluate many · P6 refuse, report or proceed, never degrade silently · P7 metadata travels with the data ·
**P8 the physics core returns true dynamics; sensors are adapters** (a_ng, τ, ω, ω̇ → observable) ·
**P9 estimation and integration are not our business** (lsqvce, gjToolbox) ·
**P10 one raise point** (`fail`), one schema per artefact (`eventTable`, `provenance`, `functionlist.csv`).

**P-11 (owner ruling 2026-09-08): two evaluation paths, both first class.** The vectorised MATLAB path (`Model.accel`, `series`, the `+kernel` units over K epochs) and the pointwise compiled path (`kernel.singleEpoch` over `mex/ngfKernel.c`, one epoch) are kept side by side. Neither is a fidelity tier of the other and neither is justified by speed; they serve an arc and a right-hand side respectively. The price of keeping both is that they must agree, and that agreement is tested at the level a caller sees, not only unit by unit.

### D.2 Layers (a layer calls only lower layers plus the four dependencies)

```
L8 app      calibrate (lsqvce bridge), forceBudget, densityInversion
L7 qc       residuals, hampel, detectImpulses, detectJumps, eventTable, screen
L6 sensor   Sensor (abstract), Electrostatic, AtomInterferometer, DragFree
L5 Model    accel | torque | series | partials | rhs | budget   <- single point of access
L4 kernel   gsiCoefficients, accommodation, dragPanelForce, radiationPanelForce,
            thermalRecoil, assemble, bodyTorques, partials
L3 geom     MacroModel, graceFO, grace, boxWing, selfShadow, thermalState, AttitudeProvider, checkAttitudeSense
L2 env      density, wind, eclipseFactor, earthRadiance, magneticField, radiationSources
L1 io       readACC1B readSCA1B readGNV1B readTHR1B readMAS1B readMAG1B, massModel,
            spaceWeather, importCERES, igrfCoefficients
L0 core     const, toolboxVersion, dependencyReport, fail, provenance, Options
```

### D.3 Workflow (GRACE-FO calibration, PDF §9)

```
L1B files ──io──▶ epoch(GPS), r,v (GCRS via RS), q(BODY→GCRS), a_meas(SRF), thr events, mass
                       │
Options ──▶ Model(macro, options, attitude, mass)
                       │  env: density+wind → v_rel; Sun+eclipse; Earth radiance → source list; IGRF
                       │  geom: normals in body frame, self-shadow fractions, panel temperatures
                       │  kernel: F_k per panel → a_ng, tau, partials  (series: chunked, parfor-optional, bit-identical)
                       ▼
sensor.Electrostatic.observe(dyn) ──▶ y_model (CoM coupling, CRN filter, frame)
                       │
qc.screen(y_meas, y_model, thr, eclipse, mass) ──▶ mask, eventTable (outlier | impulse | jump)
                       │
app.calibrate ──▶ lsqvce.fit(layout = [bias spline | 3x3 scale | jump columns], whitening from sensor.noiseModel, snoop) ──▶ Results + provenance
```
Density inversion (Doornbos) is the same pipeline with the drag scale as the unknown.

### D.4 Data types

| Concept | Type | Owner |
|---|---|---|
| epoch | referenceSystem epoch struct | RS |
| state r, v | 3×K double, GCRS, m, m/s | NGF |
| attitude | 4×K scalar-first, BODY→GCRS | NGF stores, RS operates |
| macro model | `geom.MacroModel` (validated once; packed arrays A, N, C, optVis, optIR, alpha) | NGF |
| source list | struct `.dir` 3×M body, `.P` 1×M Pa, `.band` 1×M ∈{vis,ir}, `.nu` shadow, `.meta` | NGF |
| dynamics record | struct `.a .tau .omega .omegaDot .epoch .frame .provenance` | NGF (Model → Sensor) |
| observable | struct `.value` M×K `.epoch .frame .mask` | NGF (Sensor → qc/app) |
| event table | `qc.eventTable` schema | NGF |
| series | `timetable` + provenance | NGF |

### D.5 Fidelity ladder (`Options` presets; extended models are opt-ins that no preset enables)

See `Options` help. `fast`: exponential / constant Cd / conical / no ERP / no thermal / no self-shadow; `standard`: MSIS00 / Sentman α=0.85 / conical / Knocke 12×24 / steady-state thermal / LUT; `precise`: MSIS 2.1 / Sentman+SESAM / SOLAARS-CF / CERES 24×48 / steady-state / raycast / DE440 Sun.

### D.6 Performance strategy

Vectorise over K epochs; environment in bulk (density backend accepts K points; MEX for the MSIS inner loop is the sanctioned C use); `pagemtimes` for attitude; the panel loop runs over N≈8 with K-vectorised arithmetic; self-shadow LUT replaces per-epoch ray casting; `series` splits K into chunks of `Options.series.chunk` and runs them serially or under `parfor` (PCT feature-detected) with results required bit-identical either way (test category `vectorization` + a `parallel` tag). The single-epoch C/MEX (LP-6b, delivered 2026-09-08) ports the per-epoch kernel — GSI, the gas panel loop, the radiation panel loop, the thermal recoil and `assemble` — and links nothing (see B.1: `rsApiRotation` is ITRS→GCRS and the kernel is body frame). **The justification written here until 2026-09-08 was wrong and is corrected:** it said MATLAB call overhead of 5–15 µs dominates `Model.rhs()`, but LP-6's own measurement (O-29) puts `Model.accel` at 88 ms per epoch in the standard preset and 2 ms in the fast one, so 15 µs is 0.02 % and 0.75 % — never dominant, and the cost that *is* dominant sits in `geom.selfShadow`, which this port does not touch. **(The 88 ms is superseded: WP-C3 took the standard preset to 0.798 ms per epoch at K=100, and re-measurement on 2026-09-18 found most of the 88 ms to have been per-call overhead rather than shadowing — O-29, now closed. The conclusion drawn here survives and in fact hardens: the fixed cost of one `Model.accel` call is about 15 ms in the standard preset, so D.6's 5–15 µs was wrong by three orders of magnitude in the direction that mattered.)** **Owner ruling of 2026-09-08 (P-11): the two paths are both first class.** The vectorised MATLAB path answers "what is the acceleration over this arc" and the pointwise C path answers "what is the acceleration at this one epoch, called ten thousand times"; neither is a fallback for the other and neither needs to beat the other to justify itself. The measured ratio (0.216, `tests/last_run.txt`) is therefore informational — O-31 is informational with it — and the obligation the ruling creates is agreement, not speed: `tests/tNgfMex/singleEpochEndToEndPrecision` composes the environment as `Model.evaluate` composes it, hands one epoch to the compiled kernel and requires the acceleration back to be `Model.accel`'s. The port is also a third independent implementation of the frozen kernel against the authority the mirror already is (the reason WP-33 was frozen first), and a C-callable kernel for a right-hand side that is not in MATLAB. `gpuArray` is not used (N=8 panels give nothing to a GPU; ray tracing is the only candidate and is pre-computed into the LUT).

### D.7 Environment constraints

C1 no MATLAB Coder · C2 base MATLAB + PCT only (Signal Processing and Statistics toolboxes are installed on the owner's host but **must not** be on any code path — `runNgfTests` will gain a static scan in LP-R) · C3 hand-written C only in LP-6b and the density inner loop · C4/C5 dependencies never copied · C6 references: Wöske thesis [W], HPS/XHPS [H], JPL D-56935 [J], Behzadpour 2021 [B], Siemes 2023 [S], Vielberg & Kusche 2020 [V].

---

## Part E — Standards

### E.1 Standing instruction (paste at the top of every session)

```text
=== NonGravForces STANDING INSTRUCTION ===
Follow HANDOVER.md Rev 2.0: dependency contract (Part B), principles P1-P10 (D.1),
error handling through NonGravForces.fail with the mnemonic taxonomy and a
' -> ' correction hint (E.2), arguments blocks (never inputParser), the help
template (E.3), seven test categories per public unit (E.4), no printing from
library code, full vectorisation over epochs, no aliases or wrappers for
dependency functions.
NonGravForces computes forces and torques and nothing else. referenceSystem,
gracefoReading, gjToolbox and lsqvce are called by qualified name and never
copied, vendored, forked or wrapped. Base MATLAB + Parallel Computing Toolbox
only; no MATLAB Coder; hand-written C only where C3 allows.
Deliver the function files, the test methods, the functionlist.csv row updates,
and a wrap-up with: what changed and what it was measured against; the
PREDICTED TEST COUNT; every assumption; everything not done and why.
Anything not exercised by a currently green test is a guess until the source
has been opened. If it cannot be opened, say so.
=== END STANDING INSTRUCTION ===
```

### E.2 Error handling — the design (PDF §10)

1. **One raise point.** Every error is built by `NonGravForces.fail(unit, mnemonic, message, args…)`; identifier `NonGravForces:<unit>:<mnemonic>`; the message ends with ` -> <what the caller should change>`; `fail` itself refuses a mnemonic outside the taxonomy and a message without the hint (tested).
2. **Taxonomy** (first token of the mnemonic): `input*` contract · `dep*` dependency · `data*` data file · `numeric*` NaN/Inf/singular/tolerance · `unsupported*` refused by design · `notImplemented`.
3. **Contracts first.** `arguments` blocks on every public unit (shape, class, range); semantic checks that need several inputs go through `fail` immediately after; nothing is validated twice (kernels trust packed geometry, `MacroModel` validates once).
4. **Never degrade silently.** No fallback between fidelity tiers; unknown option fields error; a missing optional dependency errors at the point of use with `depMissing` naming the setup function; `dependencyReport` warns but never refuses.
5. **Foreign errors propagate unwrapped** (`referenceSystem:*`, `gracefoReading:*`, `lsqvce:*`); NonGravForces never re-badges them.
6. **Diagnostics travel.** Every result has `info.provenance` (`provenance()`: toolbox and dependency versions, option hash, files, warnings, UTC time). Warnings issued during an evaluation are collected into `provenance.warnings`, so a post-hoc audit sees them without a console.
7. **NaN policy is explicit per unit** in the help text (e.g. `hampel`: NaN never flagged, never poisons neighbours; `density`: NaN in, `numericNaN` error, because a NaN density silently zeroes drag).
8. **Numerical guards** report, not hide: `numericSingular` carries the condition number; `numericTolerance` carries the residual; `numericOverflow` names the epoch index.
9. **Tests assert identifiers, not messages** (`verifyErrorId`), and pair every refusal with the case that must succeed (T6).

### E.3 Code rules and help template

R1 one authority per fact · R2 metadata travels with the data · R3 measured, never declared (versions read from the dependency, catalogue read from `functionlist.csv`) · R4 refuse/report/proceed · R5 do not overrule the user · R6 a tolerance is a claim about which error dominates.

Style: `camelCase` functions, `PascalCase` classes, one public unit per file, `arguments` on every public unit, SI with unit-suffixed names at interfaces (`_m`, `_K`, `_W`), radians internally, 3×K column layout, no globals, no `addpath` in package code, no printing from library code.

Help template (mandatory; enforced by `tools/scanHelp.m` in LP-R): H1 · syntax · description · INPUT ARGUMENTS (Name Size Type Units Description) · NAME-VALUE ARGUMENTS with defaults · OUTPUT ARGUMENTS · ACCURACY (every approximation named and quantified) · ERRORS (identifiers) · REFERENCES with equation numbers · EXAMPLE · See also · KEYWORDS · Created date and authoring model.

### E.4 Test rules

T1 seven categories per public unit, counted by `auditTestCategories`: `contract`, `reference`, `precision`, `vectorization`, `robustness`, `metamorphic`, `speed`. T2 exemptions in `testCategoryExemptions.m`, one reason each; `contract` never exempt. T3 a passing bound leaves its number (`verifyBoundAndRecord` → `last_run.txt`). T4 sweep before diagnosing. T5 assert ratios, not host-dependent absolutes (both are recorded). T6 informative failures, refusal paired with success. T7 default paths tested. **T8 metamorphic relations are written before the oracle** (rotation covariance, affine/shift/reversal invariance, mass and density scaling, ν∈[0,1] monotone across penumbra, series chunking invariance). **T9 an oracle is independent code** (`tools/pyxcheck/*.py`, published tables, real data), never a copy of the unit under test.

### E.5 Process rules

P-1 ☑ only after a read run report · P-2 predict then reconcile the count · P-3 freeze scope, new ideas to Part H · P-4 decisions dated · P-5 mechanisms over reminders · P-6 one substantial change per session · P-7 version rule (re-based 2026-09-06, Part I): the numeric fields of `toolboxVersion` are the authority and `semantic` is composed from them; minor = level delivered, patch = HANDOVER amendment without code change · **P-8 one zip per session, nothing beside it** · **P-9 the run log is the diary; this document holds decisions only**.

---

## Part F — The function set (authority: `doc/functionlist.csv`; rendered in `doc/html/FunctionList.html`)

| Package | Unit | Signature (frozen where ☑, proposed where ○) | WP |
|---|---|---|---|
| root | `Options(preset, name, value, …)`, `.set`, `.toStruct` | ☑ | 00 |
| root | `Model(macro, options, 'Attitude', att, 'Mass', m)` | ☑ constructor; `[a,info]=accel(ep,r,v,q)`, `[tau,info]=torque(…)`, `T=series(…)`, `[J,info]=partials(…,'WithRespectTo',[…])`, `f=rhs()`, `B=budget(…)` ○ | 14–17, 30 |
| root | `const`, `toolboxVersion`, `dependencyReport`, `fail`, `provenance`, `buildMex('Check')` (unit 52, 2026-09-07) | ☑ | 00/01, 05b |
| io | ☑ (frozen 2026-09-06) `[ep, a_srf 3×K, aux] = readACC1B(files, 'TimeRange', 'QualityBitMask')`; `[ep, q_body2gcrs 4×K, aux] = readSCA1B(files, …)` (conjugated once via `referenceSystem.invertQuaternion`, canonical q0 ≥ 0); `[ep, r, v, aux] = readGNV1B(files, 'Eop', eop|file, 'Frame', "GCRS"|"ITRS", …)` (rotated via `rotationItrsToGcrs` + `applyRotation`; `'Eop'` mandatory for GCRS); `[events, aux] = readTHR1B(files, 'TimeRange')` (`qc.eventTable` rows, class impulse, source THR1B, axis all; per-thruster on-times in `aux`); `[ep, tank, aux] = readMAS1B(files, 'TimeRange')` (gas masses, total refused ≥ 100 kg); `mass = massModel(sat, 'Mode', "constant"|"mas1b", 'Epoch', 'GasMass_kg', 'Anchor_kg')` with `mass.at(ep)`; `sw = spaceWeather(files, 'Format', "auto"|"celestrak"|"jb2008")` with `sw.query(ep)`. ○ `readMAG1B → B_meas, dipoles` (O-8); `importCERES(nc, 'Product', "SYN1deg"|"EBAF")`; `igrfCoefficients(file)` | ☑ 7 / ○ 3 | 18, 04, 08, 23 |
| env | ☑ (frozen 2026-09-07) `[rho, species 8×K, T 2×K, info] = density(model, ep, lla 3×K, sw, 'Drivers', 'AnomalousOxygen', 'Exponential', 'LocalSolarTime_h')` (exponential, msis00 via MEX; msis21/dtm2020/jb2008 refuse `notImplemented`); `[nu, info] = eclipseFactor(r, r_sun, 'Model', "cylindrical"|"conical"|"solaarsCF", 'Moon')` (solaarsCF delivered 2026-09-07 with 'EclipticNormal'; O-13 closed); `tiles = earthRadiance(backend, ep, grid, r)` (pointSource, knocke cap rings; ceres refuses); `[src, info] = radiationSources(ep, r, q_body2gcrs, options, 'Ephemeris', 'Sun', 'Moon')` with `src.dir` 3×M×K body (satellite → source), `.P` M×K Pa, `.band`, `.nu`. ○ `wind(model, ep, r, v, lla, sw)`; `magneticField(ep, r_itrs) → B` | ☑ 4 / ○ 2 | 05/25, 07, 08/27, 23 |
| geom | ☑ (frozen 2026-09-07) `MacroModel` (+ optional `vertices_m`, `heatCapacity_J_m2K`; `.V`, `.hasVertices`, `.heatCapacity`), `graceFO` (polygons, polygon centroids, Wöske capacities), ☑ `mm = grace("A"|"B", 'SolarArrayOperating')` (N-07, 2026-09-22), `mm = boxWing([Lx Ly Lz], optics, 'Array', 'Name', 'Alpha', 'HeatCapacity_J_m2K')`, `[f NxK, info] = selfShadow(macro, dir_body 3xK, 'Method', "raycast"|"lut", 'Samples', 'LutResolution', 'Lut')` (dir from the body TO the source; f = 0 on unlit plates), `[T NxK, info] = thermalState(macro, src, 'Model', "steadyState"|"transient", 'Exposed', 'InternalPower_W', 'InternalShare', 'Epoch', 'T0_K')`, `ap = AttitudeProvider(epochNodes, qNodes, 'Method', "slerp"|"squad")` with `[q, info] = ap.at(epoch)` and the static `AttitudeProvider.nominal(r, v, xAlongVelocity)`, `[ok, report] = checkAttitudeSense(q, r, v, 'NadirTolerance_deg')` | ☑ 7 | 09–11, 28 |
| kernel (LP-6b) | ◐ `[F 3×3, tau 3×3, info] = singleEpoch(macro, gas, src, Tpanel_K, 'Model', 'Alpha', 'Tw_K', 'Cd', 'SpecularFraction', 'CdScale', 'CrScale', 'CentreOfMass')` — one epoch through the compiled kernel; columns of `F` and `tau` are drag \| radiation \| thermal, `info.Fk` is 3×N×3. Refuses `depMissing` until `NonGravForces.buildMex` has run | ◐ | 17b |
| kernel | ☑ (frozen 2026-09-07) `[Cp, Ctau, info] = gsiCoefficients(model, s, theta, alpha, Tw_K, Tinf_K, 'Cd', 'SpecularFraction')` (normal and tangential per unit plate area, (1/2) rho U^2 convention; sentman, maxwell, constantCd; cll/schamberg refuse, O-17); `[alpha, info] = accommodation(model, species 8×K, T 1×K, 'Alpha')` (fixed; sesam/langmuir refuse, O-16); `[F 3×N×K, info] = dragPanelForce(macro, rho, species, T, vrel_body, 'Model', 'Alpha', 'Tw_K', 'Cd', 'SpecularFraction', 'Shadow')`; `[F, info] = radiationPanelForce(macro, src, 'Shadow' N×K×M)`; `[F, info] = thermalRecoil(macro, T_K N×K)`; `[F, tau, info] = assemble(Fk, macro, 'CentreOfMass')`; `[tau, info] = bodyTorques(macro, r_body, 'GravityGradient', 'Magnetic', 'B_body', 'Dipole_Am2', 'GM')`; `partials(…)` ○ (LP-6) | ☑ (partials ○) | 12, 13, 26, 24, 16 |
| sensor | `Sensor` ☑ (abstract: `observe`, `partials`, `noiseModel`, `numParameters`); `Electrostatic(bias, S, dr_proofMass, filter)`; `AtomInterferometer(k_eff, T, cycle, axis, 'Rotation', …)`; `DragFree(controller)` | ○ | 29 |
| qc | ☑ (frozen 2026-09-07) `hampel`, `eventTable`, `[r, info] = residuals(y_meas K×M, y_model, 'Bias', 'Scale', 'Whitening' none|difference|highpass, 'Window')`; `[events, info] = detectImpulses(r, 'Epoch' K×1 datetime UTC, 'Response', 'Priors', 'Threshold', 'PriorWindow_s')` (zero-mean matched filter); `[events, info] = detectJumps(r, 'Epoch', 'Model' level|slope, 'Penalty', 'MinSegment', 'Priors', 'PriorWindow_s')` (PELT, BIC (p+1) sigma^2 log K); `[mask, events, info] = screen(y_meas, y_model, epoch, …, 'Response', 'Thrusters', 'Priors', 'JumpMinSegment')` | ☑ | 31 |
| app | `[cal, info] = calibrate(yMeas, aModel, epoch, 'Bias', "constant"\|"linear"\|"hat", 'KnotInterval_s', 'Scale', "none"\|"diag"\|"full", 'Events', 'Weights', 'Mask', 'Robust', "huber")` — the toolbox's own separable weighted least squares, no external estimator; `forceBudget(Model, ep, r, v, q) → table`; `densityInversion(…)` | ○ | 19, 30, 32 |

Counting: 8 root + 10 io + 6 env + 7 geom + 9 kernel + 4 sensor + 6 qc + 3 app = **53 units** (`buildMex` added 2026-09-07, `kernel.singleEpoch` 2026-09-08), of which 43 ☑/◐/interface, 10 ○ stubs (the catalogue is the authority: 42 done + 1 interface + 10 stub) (the 39th stub row, `Sensor`, became the interface; LP-2 turned seven io stubs into units on 2026-09-06). Every unit has one catalogue row, one file, one `Contents.m` line, and — when ☑ — seven test categories or a reasoned exemption.

Simplifications agreed by the expert panel (Part J): one density dispatcher instead of three classes; one `eclipseFactor` with a model option instead of three functions; the source list makes SRP, albedo and IR one kernel path; sensors replace per-instrument calibration code; no compat layer; no gravity; no EOM; no integrator; no least squares.

---

## Part G — Level prompts (the plan; one level = one session; WP numbers preserved where content survived)

Each task states: what · what must be green first · deliverables · accuracy as a number with its reason · **the oracle** · re-entry criteria.

### LP-0 — Foundation ☑ (this session)
WP-00 skeleton, `const`, `Options`, `fail`, `provenance`, catalogue, stubs · WP-01 `dependencyReport` (measured) · WP-03 test infrastructure (seven categories, audit, exemptions, one runner). Oracles: CODATA 2018, IAU 2015 B3, JPL D-56935. Result: 74/74 GREEN.

### LP-1 — does not exist (time and frames are referenceSystem).

### LP-2 — Data I/O ☑ (2026-09-06; run report 116/116 GREEN, count predicted 116)
WP-18 readers over `gracefoReading` (GNV1B rotated ITRS→GCRS through `rotationItrsToGcrs`; SCA1B conjugated once; THR1B as event table; MAS1B tank masses; MAG1B) · WP-04 `spaceWeather` (CelesTrak SW-All + SET SOLFSMY/DTCFILE for JB2008; observed and previous-day F10.7; ap7 history; Kp; Dst) · `massModel`. Oracles: the owner's real 2020-01-01 GRACE-C samples; handbook mass; a Greenwich-angle test of the rotation direction. Accuracy: GNV1B round trip < 1 mm; time tag round trip exact.
**Delivered 2026-09-06**: seven units (MAG1B stays a stub, O-8), private helpers `io/private/l1bLoad`, `l1bColumns`, test class `tNgfIo`, real-data fixtures (Part I). Measured on the 2020-01-01 GRACE-C sample: GNV1B round trip 5e-9 m; 5-point finite-difference velocity vs rotated velocity 4.4e-5 m/s (4.7e2 m/s without the `Rdot` term); GCRS−ITRS longitude of record 1 vs the IAU 1982 GMST polynomial −0.27° (precession since J2000); body +Z vs nadir 0.86° conjugated, 96° not. Result: 116/116 GREEN (Tier A, host MATLAB R2026a Update 4), Code Analyzer 0 findings, `auditPackage.py` clean; first run had 6 red rows, all one cause (O-12) plus one test regexp — decomposed before any bound moved.

### LP-3 — Environment ☑ (2026-09-07; run report 142/142 GREEN, count predicted 142)
WP-05 `density` with backends exponential and msis00 (MEX of the Brodowski C, oracle = its 17-case DOCUMENTATION table, reproduced exactly) · WP-07 `eclipseFactor` (cylindrical, conical, solaarsCF; oracle Montenbruck & Gill §3.4.2 example and a numerical disc-integration in Python; C¹ test across penumbra) · WP-08 `earthRadiance` (pointSource, knocke; oracle Knocke 1988 Table 1 and Orekit's `KnockeRediffusedForceModel` numbers) + `radiationSources` (metamorphic: rotating attitude rotates the source list). Accuracy: density backend bit-exact to the C reference; ν to 1e-12; Knocke flux to 1e-6 relative.
**Scope frozen 2026-09-07:** `density` with `exponential` (caller-parameterised, sizing only) and `msis00` (MEX of the vendored Brodowski C, `gtd7d` by default, `gtd7` via `'AnomalousOxygen', false`; `msis21`, `dtm2020`, `jb2008` refuse with `notImplemented`, LP-3b) · `eclipseFactor` with `cylindrical`, `conical` and the lunar factor (`'Moon'`); `solaarsCF` refuses with `notImplemented` (O-13) · `earthRadiance` with `pointSource` (Lambert sphere, far field) and `knocke`; `ceres` refuses with `notImplemented` (WP-27) · `radiationSources` over `Options` (`rs.analytic`; `rs.de440` through `'Ephemeris'`, refusing with `depKernel` when none is given) · new public unit `buildMex` (WP-05b) · oracles: DOCUMENTATION §7 table (17 cases) for msis00; pure-Python disc integration `tools/pyxcheck/eclipse_check.py`; pure-Python cap integration `tools/pyxcheck/knocke_check.py`; Orekit numbers are NOT in hand and are not used.
**Delivered 2026-09-07 and measured** (every number recorded per run in `tests/last_run.txt`): the 17 DOCUMENTATION cases reproduce digit for digit (worst 4.8e-7 = the printed 7th digit); conical shadow vs the Python integration 5.6e-14 after the atan2 rewrite (6.4e-10 with the acos form: conditioning finding, Part I); rotation invariance 2.5e-14; Knocke cap tiles vs the 0.25° Python integration 2.3e-2 at [12 24], 5.4e-3 at [24 48], 1.4e-3 at [48 96] (first-order convergence, ratio 4.2); point source at 100 R: albedo 7e-3, longwave 3e-5; SRP at 1 au = tsi/c exactly; the standard source list for the 600-epoch real sample 0.07 s (M = 577). First run had 8 red rows: 5 test-authoring slips, one Options validator message, two unit findings (eclipse NaN not refused; acos conditioning) — decomposed, no bound moved without a measured reason.

### LP-3b — Extended environment
WP-25 `msis21`, `dtm2020`, `jb2008` backends (oracles: NRL/CCMC and CNES test outputs; SET JB2008 test cases), `wind` with HWM14 (oracle NRL test values) · WP-23 `magneticField` IGRF-14 (oracle igrf14.f test cases; coefficients from Zenodo 10.5281/zenodo.14012303) · WP-27 CERES import and ADM anisotropy (oracle: a real SYN1deg file header — **open item B7 from v0.1.1 is still open**: the netCDF variable names must be read from a real file before this WP is written).

### LP-4 — Geometry and attitude ☑ (2026-09-07; run report 169/169 GREEN, count predicted 169)
WP-09 `boxWing`, cylinder/sphere kinds (the S-band boom) · WP-10 `selfShadow` raycast + LUT with a convergence sweep over `lutResolution` (T4) · WP-11 `AttitudeProvider` over `interpolateQuaternions`, `checkAttitudeSense` by physics (nadir, line of sight) · WP-28 `thermalState` steady-state and transient (oracle: [W] Chapter 4 figures re-digitised; energy balance metamorphic: total absorbed = total emitted in steady state).
**Scope frozen 2026-09-07:** row 0: `eclipseFactor` model `solaarsCF` (closes O-13; Robertson 2015 thesis Sect. 3.7, Table 3.2, transcribed 2026-09-07, see Part I) · WP-11 `AttitudeProvider(source, …)` over `referenceSystem.interpolateQuaternions` (SCA1B nodes; `nominal` from `rotationGcrsToRsw` and `quaternionFromRotation`) with `.at(epoch)`, refusing extrapolation as the dependency does; `checkAttitudeSense(q, r, v)` by physics (body +Z within a tolerance of nadir, body +X along ±velocity) · WP-09 `MacroModel` gains optional plate `vertices` (3xV per panel, used only by self-shadowing) and `boxWing(dims, optics)`; the S-band boom stays out (its axis is in no source in hand, R3); `cylinder`/`sphere` kinds stay reserved (no kernel consumes them before LP-5) · WP-10 `selfShadow(macro, dir_body, 'Method', "raycast"|"lut")`: stratified ray cast on the plate polygons, icosphere look-up table over incidence directions, convergence sweep over `lutResolution` (T4) · WP-28 `thermalState(macro, src, 'Model', "steadyState"|"transient")`: steady state from the absorbed flux balance; transient per-panel lumped ODE (Wöske 2019 eq. 5.18) with per-panel heat capacity per area from [W] Table 6.1 layers; oracles: energy balance (steady), linearised exponential relaxation (transient), pure-Python ray cast for `selfShadow`, the real SCA1B day for attitude.
**Delivered 2026-09-07 and measured** (per run in `tests/last_run.txt`): ray cast vs the pure-Python cast 1.25e-2 at 40 samples, 3.8e-3 at 100 (sampling bound ~2/m); two-plate analytic case exact; LUT mean error 0.089 → 0.056 → 0.026 over 42/162/642 vertices, LUT query 0.19 of the ray-cast time; nominal attitude within 0.86° of the real SCA1B day, the other sense 180° off; steady-state energy balance 6e-16; transient vs linearised relaxation 1.1e-3 K; the LP-0 flank centroids were 0.39 m wrong (Part I); the first run had 5 red rows: one LP-0 test that predates polygons, three test slips (a bitwise expectation after renormalisation, a composition order, an acos-conditioned bound) and one unit finding (batch-size-dependent BLAS summation in the shadow fraction, made column-wise) — no bound moved without a measured reason; a second unit finding followed (fully covered plate 1e-16 instead of 0, made exact).

### LP-5 — Kernel and quality control ☑ (2026-09-07; run report 231/231 GREEN, count predicted 231)
**Scope frozen 2026-09-07 (Tier A, MATLAB MCP bridge over the repo; CI tier through the PR):** `[Cp, Ctau] = gsiCoefficients(model, s, theta, alpha, Tw_K, Tinf_K, 'Cd', 'SpecularFraction')` returns the NORMAL and TANGENTIAL coefficients per unit plate area in the (1/2) rho U^2 convention (sentman: incident Maxwellian flux plus diffuse re-emission at the Koppenwallner/Sentman energy-equivalent temperature; maxwell: a specular fraction that doubles the incident normal momentum and carries no shear or re-emission; constantCd: Cp = Cd*max(cos theta, 0), Ctau = 0, so the plate sum reproduces eq. 5.20 of [W] with the projected area; cll and schamberg REFUSE `unsupportedModel` — no closed form and no oracle in hand, see O-17) · `alpha = accommodation(model, species, T, 'Alpha')` (fixed; sesam and langmuir refuse `notImplemented`: the SESAM constants of Pilinski 2010/2013 are not in hand and are not guessed, O-16) · `[F, info] = dragPanelForce(macro, rho, species 8xK, T 1xK, vrel_body 3xK, 'Model', 'Alpha', 'Tw_K', 'Cd', 'SpecularFraction', 'Shadow')` → 3xNxK with per-species mixing of the coefficients (s_j from each molar mass, const.speciesMass) · `[F, info] = radiationPanelForce(macro, src, 'Shadow')` (theory eq. 6 = [W] eq. 5.4, vis and ir triplets by `src.band`) · `[F, info] = thermalRecoil(macro, T_K NxK)` ([W] eq. 5.14) · `[F, tau] = assemble(Fk, macro, 'CentreOfMass')` · `[tau, info] = bodyTorques(macro, r_body, 'GravityGradient', 'Magnetic', 'B_body', 'Dipole_Am2', 'GM')` (gravity gradient 3 GM/r^5 r x I r; magnetic m x B; both opt-in logicals, nothing enabled by default) · `[r, info] = residuals(y_meas, y_model, 'Bias', 'Scale', 'Whitening', 'Window')` (the sensor object of the Part F sketch arrives with LP-7; the preliminary calibration is passed explicitly) · `[events, info] = detectImpulses(r, 'Response', 'Epoch', 'Priors', 'Threshold', 'Window')` (matched filter; THR1B priors assigned first) · `[events, info] = detectJumps(r, 'Epoch', 'Model', "level"|"slope", 'Penalty', 'MinSegment', 'Priors', 'PriorWindow_s')` (PELT, Killick 2012, BIC penalty) · `[mask, events, info] = screen(y_meas, y_model, epoch, 'Bias', 'Scale', 'Response', 'Thrusters', 'Priors', ...)`. `partials` (WP-16) stays with LP-6 as Part G says. Oracles: `tools/pyxcheck/sentman_check.py` (numerical kinetic-theory integration, written 2026-09-07 before the unit), the hyperthermal and static-gas limits in closed form, the Lambert sphere Cr = 1 + 4/9 delta and the absorbing/specular spheres, Stefan-Boltzmann in closed form and the steady-state balance with `thermalState`, the gravity-gradient closed form for a diagonal inertia, `tools/pyxcheck/pelt_check.py` (exhaustive optimal segmentation on short planted series), synthetic pulses at known epochs with the real THR1B day as priors, and the real 2020-01-01 ACT1B sample for the screen pipeline. [S] Table 3 and [V] monthly means need the whole model (LP-6): O-6 stays open.
**Delivered 2026-09-07 and measured** (per run in `tests/last_run.txt`): Sentman Cp/Ctau vs the kinetic-theory quadrature 4e-14 / 1.2e-14 over 105 cases, the drag and lift projections 4e-14 / 1e-14, the specular case 5e-14; Cp gamma + Ctau l equals Doornbos' C_D to 1e-15; hyperthermal limit 2 gamma^2 up to the re-emission term 9.71e-5 at s = 1e4; GRACE-FO Sentman Cd on the projected area 3.59 (published band 3.2-3.8); constantCd plate sum vs eq. 5.20 to 1e-16; Lambert/absorbing/specular spheres of 4000 plates Cr within 4e-5 of 1 + 4/9, 1, 1; steady-state recoil = (2/3) absorbed/c to 1e-12 (was 7e-7 before O-2); gravity gradient vs Wertz 17-30 to 1e-14; PELT equals the exhaustive optimum on 6 (model, penalty) cases bitwise; planted 4-sigma step to +-1 sample and 3 sigma; a 0.1 sigma/sample kink to 4 samples; planted CRN pulses to +-1 sample and 3 sigma_a; the real THR1B day assigned as priors; the real 600-s ACT1B sample screened (7 impulse candidates recorded, no oracle). The first kernel run had 7 red rows: two unit precision findings (1 + erf cancelling below zero on leeward plates -> erfc; the alpha bracket at alpha = 1 -> rewritten), one unit contract defect (thermalRecoil refused the 0 K of unlit steady-state plates), four test slips (an edge-on plate does feel 1/(2 s^2), two hand numbers wrong, an expansion-compatible shape). The qc run had 4 red rows: BIC parameter count (unit), Hampel-flattened pulse peaks before the matched filter (unit, screen order), two guessed bounds restated from measurements. PELT bug found by the oracle before the tests: candidates pruned by omission. No bound moved without a measured reason.
**Finding made while reading (2026-09-07):** the reference document (`theory.txt` §2.3 eq. 4-5) labels Doornbos' flat-plate DRAG and LIFT coefficients as `Cp` and `Ctau` and then uses them in eq. 1 along the normal and tangential directions; in the hyperthermal limit that form gives a normal coefficient 2 gamma instead of 2 gamma^2. The kernel derives the normal and tangential coefficients from the momentum integrals (Cp_inc = gamma P/sqrt(pi) + (gamma^2 + 1/(2 s^2)) Z, Ctau_inc = l (P/sqrt(pi) + gamma Z)); Cp gamma + Ctau l reproduces Doornbos' C_D and Cp l - Ctau gamma his C_L, which the tests assert. The document is to be corrected at LP-R (WP-21).
WP-13 `gsiCoefficients` Sentman exact (oracle: Python kinetic-theory integration in `tools/pyxcheck/sentman_check.py`; the v1 ×8 defect as a regression test; flat-plate limits s→∞) · WP-26 CLL, Maxwell, Schamberg, `accommodation` SESAM (oracle Pilinski 2013 Fig. 3 values) · WP-12 `radiationPanelForce`, `thermalRecoil` (oracle: [V] monthly SRP/ERP means; [S] Table 3 RMS per axis within 20 %) · WP-24 `bodyTorques` (oracle: gravity-gradient closed form for a diagonal inertia; magnetic torque m×B) · WP-31 `qc`: `residuals`, `detectImpulses` (oracle: synthetic CRN-shaped pulses at known epochs, THR1B priors from the real file), `detectJumps` (oracle: planted steps, PELT vs brute-force segmentation on short series), `screen` (metamorphic: time shift, affine scaling, reversal). Accuracy: kernel partials vs complex step 1e-8; impulse epoch to ±1 sample; jump amplitude to 3σ.

### LP-5x — Cross-validation tooling ☑ (2026-09-07; run report 238/238 GREEN, count predicted 238; WP-34 not delivered, O-19)
**Scope frozen 2026-09-07 (Tier A, MATLAB MCP bridge over the repo plus `matlab -batch ciGate`; CI tier through the PR).** WP-33 is delivered as a **standard-library** mirror, not a numpy one: the host has Python 3.14.7 and no numpy (measured 2026-09-07, `ModuleNotFoundError`), and the ubuntu audit job installs nothing, so a numpy mirror could be regenerated on neither tier. Deliverables: (a) `tools/pyxcheck/ngfkernel/` — flat stdlib modules `gsi.py`, `drag.py`, `radiation.py`, `thermal.py`, `assemble.py`, `torques.py` mirroring the frozen kernel (`gsiCoefficients` sentman/maxwell/constantCd, `dragPanelForce` with per-species mixing, `radiationPanelForce`, `thermalRecoil`, `assemble`, `bodyTorques`; `accommodation` fixed is a pass-through and is mirrored inside the driver). Written from the equations by agents instructed to open none of `+NonGravForces/+kernel/`, `tools/pyxcheck/sentman_check.py`, `tests/`, `HANDOVER.md` or `CLAUDE.md`. **Measured, not assumed (correction to this row, 2026-09-07): the exclusion of `CLAUDE.md` was NOT achievable** — the harness injects it into every subagent automatically, and several agents disclosed that unprompted. Its one implementation-adjacent line is the Sentman pitfall, which the derivations CONTRADICTED (see Part I), so the leak pushed the wrong way. Two of the fourteen agents also set a self-reported contamination flag while naming only permitted sources in their prose, so that flag is not evidence either way. Independence therefore rests on the arithmetic, not on trust: the discriminating edge-on limit, the disagreements with the reference document, and the parity numbers below (T9); the hardest quantity (the gas-surface coefficients) is derived three times by three routes and the routes are reconciled numerically before the module is accepted. (b) `tools/pyxcheck/kernel_mirror_check.py` → `kernel_mirror_oracle.csv` over a committed case list, with `kernel_mirror_plates.csv` (the GRACE-FO and cube plate tables exported once from `geom`) as INPUT data: `geom` is not the unit under test, the ruling mirrors the kernel only. (c) `tests/tNgfMirror.m`: MATLAB versus mirror on every case, per unit. (d) `tools/pyxcheck/hampel_check.py` converted to the standard library, reproducing the committed `hampel_oracle.csv` bit for bit (found 2026-09-07: it imports numpy, so the one committed oracle in the tree could not be regenerated on either tier — an oracle that cannot be re-run is a claim, not an instrument). **WP-34 is not delivered: GROOPS is not installed on the host (searched 2026-09-07) and is the owner's to supply (O-19).** After this level the mirror is FROZEN: a regression oracle, not a maintained second product.
**Delivered 2026-09-07 and measured** (per run in `tests/last_run.txt`). Parity of the MATLAB kernel against the mirror, relative to the scale of each case: `gsiCoefficients` cp 8.7e-14 over 617 of 760 cases and ctau 6.0e-15 over 468 (the remaining cases are leeward values below 1e-12, where both sides are required to be negligible and are, 3.7e-13); `dragPanelForce` 5.9e-16 over 80 cases; `radiationPanelForce` 2.3e-16 over 16 cases and its absorbed flux 4.8e-16; `thermalRecoil` 2.1e-16; `assemble` and `bodyTorques` exactly 0. The exported plate fixture still matches live geom to 4.4e-16. `hampel_check.py` reproduces its committed flags bit for bit with the standard library (400 rows, 5 flagged). **The headline is a null result**: three independent derivations of the gas-surface coefficients reached algebraically identical forms, and NONE reproduced the drag/lift mislabel of the reference document — the trap that LP-5 had found by a different route. The discriminating limit is edge-on incidence, where the true normal coefficient is 1/(2 s^2) and the tangential one l/(s sqrt(pi)) with no 1/s^2 term at all; both sides reproduce both to 5e-16 at s = 0.5, 2 and 8. The planned completeness-critic stage did not run (its agent died on a session limit); the audit was done by hand instead and is recorded in Part I.
**WP-33 `tools/pyxcheck/ngfkernel/`** — a numpy mirror of the FROZEN KERNEL ONLY: `gsiCoefficients`, `dragPanelForce`, `radiationPanelForce`, `thermalRecoil`, `assemble`. Written from the equations, not transcribed from the MATLAB, or it inherits the author's misunderstandings and stops being an oracle (T9). No mirror of io, env, geom, Model, sensor, qc. **Frozen when LP-5 is GREEN**: from then on a regression oracle, not a maintained second product. Rationale: LP-6b ports this same kernel to C, so both ports are then checked against one authority.

**WP-34 GROOPS bridge** — GROOPS is used as a BLACK BOX (GPLv3: no source is transcribed into this toolbox in any language). `SimulateAccelerometer` with `miscAccelerations:atmosphericDrag` (Sentman 1961 in the Moe & Moe 2005 analytical form, `considerTemperature`, `considerWind` HWM14) and `miscAccelerations:radiationPressure` (SRP + albedo/longwave grids + TRP) over the real 2020-01-01 GRACE-C day. Outputs are frozen as CSV under `tests/fixtures/groops/` and committed; neither the CI nor the owner's host needs a GROOPS installation. **Independence caveat:** GROOPS cites Wöske et al. (2019) for radiation and TRP — the same source as reference [W] — so for SRP/ERP/TRP this is an implementation comparison, not an independent oracle; the drag path (Moe & Moe analytical vs our exact Sentman) is a genuine cross-check.

Side effect: `gracefoReading.writeGroopsInstrument` / `readGroopsInstrument` document that their header conventions were never verified against a real GROOPS installation. WP-34 verifies them and reports back to that toolbox.

### LP-6 — Model and partials ☑ (2026-09-08; run report 262/262 GREEN, count reconciled in Part I)
**Scope frozen 2026-09-08 (Tier A: MATLAB MCP bridge for single classes, `matlab -batch ciGate` for the suite; CI tier through the PR).**
Delivered here: `env.wind` (WP-25, corotation only), `kernel.partials` (WP-16) and the six `Model` evaluation methods (WP-14, WP-15, WP-17).
`app.forceBudget` (WP-30) moves to LP-7 — see the decision log; it is the application layer and a thin reduction over `Model.budget`, which lands here.

**Interfaces frozen.**
`[vrel, info] = env.wind(model, epoch, r_gcrs, v_gcrs, 'Eop', eop, 'Wind_m_s', w)` → 3×K GCRS velocity of the satellite relative to the air.
`none` returns v; `corotation` returns `v - Rdot*(R'*r)` with R and Rdot from `referenceSystem.rotationItrsToGcrs`, so the co-rotation term is
the dependency's, never a hand-written ω×r (theory §2.2); `hwm14` refuses `notImplemented` (O-27). Measured on the 2020-01-01 arc: 255.8 m/s at −59.5° latitude.
`[dF, names, info] = kernel.partials(what, macro, ...)` → analytic per-panel derivatives of the drag and radiation panel forces with respect to the
parameters that enter their closed forms.
`Model(macro, options, 'Attitude', 'Mass', 'Eop', 'SpaceWeather', 'Ephemeris')` — the constructor gains the last three: `accel` needs the geodetic
height (hence the ITRS rotation, hence EOP) and the density drivers, and a Model that lacks them must refuse by name rather than assume a default.
`[a, info] = accel(ep, r, v, q, 'Frame', "body"|"gcrs")` (body by default, principle P4), `[tau, info] = torque(...)`,
`T = series(ep, r, v, q, 'Quantities', ...)` (chunked by `options.series.chunk`, `parfor` optional, bit-identical to the unchunked result),
`[J, names, info] = partials(ep, r, v, q, 'WithRespectTo', ["r" "v" "cdScale" "crScale" "alpha"])`, `f = rhs()`, `B = budget(...)`.

**The complex-step verifier the line above asks for is impossible on this path, measured 2026-09-08**: MATLAB's `erf` and `erfc` reject complex
input (`MATLAB:erf:notFullReal`), and `abs`/`vecnorm` discard the imaginary part, so no complex step can pass through Sentman's coefficients or any
norm in the chain. The verifier is Richardson-extrapolated central differences instead, and the accuracy claim is stated against that.
`partials` is analytic and exact where a parameter enters linearly (`cdScale`, `crScale`), analytic through the kernel closed form for `alpha`, and
central-difference for `r` and `v`, whose chain runs through the MSIS MEX; `Options.partials.method` selects, and `complexStep` refuses with the
measurement above as its reason.
`rhs()` was delivered as code and stayed ◐ against its consumer for two levels; in LP-C (WP-C9) it is **removed** and `app.DynamicsSource` serves gjToolbox's frozen interface instead (O-4, issue #35). `Model`'s methods are now `accel | torque | series | partials | budget`.

**Oracles (none built here):** the real 2020-01-01 GRACE-C day — ACT1B along-track acceleration against `Model.accel`, the drag direction against
−v_rel, and the eclipse against the shadow function; published magnitudes for the thermal recoil ([S] Siemes 2023, 6.4 nm/s² cross-track RMS;
the dry run of 2026-09-08 gives 6.1 nm/s²); Richardson-extrapolated central differences for every analytic partial; the unchunked series for the
chunked one, bitwise; `budget` summing to `accel` bitwise; the independent Python kernel mirror already frozen at LP-5x underneath it all.
**Dry run 2026-09-08 before freezing** (5 epochs of the real day at 523 km, ρ = 5.06e-14 kg/m³, 577 radiation sources, ν = 1):
drag 8.5e-9 m/s² (7.0e-9 along −v_rel), radiation 2.6e-8, thermal 6.1e-9, total 2.7e-8 m/s². The chain composes end to end today.
**Delivered 2026-09-08 and measured** (per run in `tests/last_run.txt`). Against the real 2020-01-01 GRACE-C day, standard preset: the modelled acceleration correlates with the measured ACT1B at **0.980 cross-track and 0.994 radial** over 200 s; along-track it is 0.190 and is RECORDED, not asserted, because the ten-minute excerpt gives a measured 3.2e-7 RMS against a modelled 9.4e-9, dominated by the instrument's own bias drift. Thermal recoil median 7.9e-9 m/s² against the 6.4 nm/s² cross-track RMS Siemes et al. (2023) publish for GRACE. The analytic parameter partials agree with central differences in the option to 4.6e-10 (cdScale), 1.2e-10 (crScale) and 5.5e-10 (alpha); the state partials predict a 37 m step to 3.3e-8 relative. `kernel.partials` agrees with Richardson-extrapolated central differences to 1e-12 (alpha), 7e-11 (speed), 2e-11 (densityScale), 8e-11 (pressureScale) and 2e-11 (temperature). accel minus its three effects, budget against accel, and the rhs handle against accel are all bitwise zero; a series in chunks of two equals the unchunked one bitwise; body and GCRS magnitudes agree to 3e-24. `env.wind`: the co-rotation follows omega*sqrt(x²+y²) to 2e-6 across an arc spanning 256 to 467 m/s. Cost measured: 88 ms per epoch in the standard preset (0.5 s once at construction for the 642-vertex table), 2 ms per epoch in the fast preset (O-27). **Both figures are of 2026-09-08 and are superseded by WP-C3 (0.798 ms and 0.076 ms at K=100, measured 2026-09-18); the 88 ms was also a K=1 measurement, where a ~15 ms fixed per-call cost dominates — O-29.** **Finding made by running the P-11 end-to-end test (2026-09-08):** `kernel.singleEpoch` screened `gas.T_K` for finiteness whatever the gas-surface model, and `env.density("exponential")` returns species AND temperature as NaN **on purpose**, because it models neither. That is exactly the pairing the fast preset uses (exponential density with constantCd), so the compiled path refused the fast preset's own environment and nobody had noticed: every earlier test built its gas by hand. An input a model does not read is now not validated, the refusal is kept for the models that do read it, and both halves are asserted. This is the argument for composing a test's inputs from the real units rather than from literals. **A second bound in the same test was mine and wrong the same way as the LP-6b speed bound:** it claimed the relative velocity would match the Model's bitwise, and the runner measured 9.095e-13 m/s. `Model.evaluate` rotates with `pagemtimes` over K pages, the test with one matrix-vector product, and a BLAS product sums in a batch-size-dependent order — a pitfall this document has carried since LP-6 and the test author did not apply. 9.1e-13 on 7600 m/s is 1.2e-16 relative; the claim is now relative agreement with the routine difference named.

**Findings made by running:** the drag is NOT antiparallel to the stream — Sentman's shear is genuine lift, 0.34 deg here, while constantCd gives 1e-16 by construction, so the test asserts the second and measures the first; and `geom.selfShadow`'s look-up interpolation could select a face the ray is nearly parallel to and return NaN, which no test could see (Part I).
WP-14 `Model.accel/torque` (integration requirement I2 from v0.1.1 kept: the kernel returns per-panel forces, never sums internally) · WP-15 `series` chunked/parallel with bit-identity test · WP-16 `partials` (r, v, Cd, Cr, α) · WP-17 `rhs()` and `budget` hook for gjToolbox · ~~WP-30 `forceBudget`~~ (moved to LP-7, 2026-09-08). Oracles: [S] and [V] magnitudes; ACT1B along-track drag at solar minimum ~5e-8 m/s².

### LP-6b — Single-epoch C/MEX ☑ (2026-09-08; written on Tier B, verified through CI on the self-hosted runner: 270/270 GREEN after one guessed bound was decomposed)
WP-17b transcribes the frozen kernel to C for ONE epoch: `gsiCoefficients`, the gas panel loop, the radiation panel loop, `thermalRecoil` and `assemble`, plus the composite `ngfSingleEpoch` that applies `cdScale`/`crScale` where `Model.evaluate` applies them. `mex/ngfKernel.h`/`.c` is C99 with no MATLAB headers, no allocation, no globals and no printing; `mex/ngfKernelGateway.c` marshals and validates shapes; every constant travels in from MATLAB (`const`, `iersConstants`), because a constant written twice is a constant that will disagree once. MATLAB entry point `kernel.singleEpoch`, binary in `+kernel/private`, built by the second `buildMex` target.

**Two scope rulings, taken in session and APPROVED by the owner on 2026-09-08.** (i) `thermalRecoil` is ported although the prompt named only GSI + radiation + assemble: the frozen mirror covers it, it is six lines of C, and without it the MEX cannot return the acceleration `Model.accel` returns, which is the whole point. (ii) `rsApiRotation` is NOT linked — see B.1; it is ITRS→GCRS and the kernel is body frame throughout, so there was nothing for it to do. The line in the prompt survived from a planning revision and is now corrected in B.1 and D.6.

**Measured 2026-09-08 without MATLAB** (`bash mex/test/runCChecks.sh`, gcc 13.3, `-Wall -Wextra -Werror`), the C against the frozen Python mirror with the bounds `tests/tNgfMirror.m` applies to the MATLAB: `gsiCoefficients` 8.732e-14 relative on the 617 coefficients above 1e-12 and 6.034e-15 on the 468 tangential ones; `dragPanelForce` 4.328e-16; `radiationPanelForce` 1.990e-16 and its absorbed flux 2.701e-16; `thermalRecoil` 1.388e-16; `assemble` force and torque exactly 0; the composite equals the four units called separately bitwise; 11 refusal guards provoked and all behaved. **The two gsi figures are the MATLAB's own numbers against that mirror to five significant digits** (8.7315e-14 and 6.0339e-15 in `tests/last_run.txt`), i.e. the residual is the mirror's better-conditioned leeward form (O-20) and not a difference between the two ports; transitively |C − MATLAB| ≲ 1e-15 on the forces. That is a BOUND, not a measurement: the measurement is `tests/tNgfMex.m` and it needs MATLAB. The harness was fault-injected (a 1e-12 relative error in the Lambert factor) and turned red, so it has teeth.

**Verified in MATLAB 2026-09-08 through the CI `matlab` job** (the MCP bridge never answered all session; the self-hosted runner did, which is how a Tier B session got a run report at all — see the note under A.3). Predicted 270, ran 270: **269 passed, 1 failed** on the first run, **270/270 GREEN, 0 incomplete, 0 audit shortfalls, gate TRUE** on the second once the guessed bound below was decomposed (R2026a Update 5, 218.9 s, log in `tests/last_run.txt`). Both gateways built with MinGW64 14.2.0. THE 1e-14 CLAIM IS NOW MEASURED, not bounded: force **4.4897e-16**, torque **1.9492e-15**, per-panel gas force **5.9779e-16** against the MATLAB kernel chain over five epochs; `gamma` and the per-species speed ratio agree **bitwise** (0.0); mixture `Cp` 3.1496e-15. Through the composite against the mirror: gas force 4.3992e-16, radiation 3.4502e-16, absorbed flux 4.7559e-16, thermal recoil 6.9807e-15.

**Speed, measured** (8 plates, 3 radiation sources, median of 5 x 200 calls on the runner): compiled **901.8 µs** per epoch against **4149.6 µs** for the MATLAB kernel chain computing the same three effects, ratio **0.217** (the green rerun measured 836.1 against 3879.0, ratio 0.216 — the ratio is stable to 0.5 % across two runs, the absolutes are not, which is the whole reason T5 asks for the ratio) — the compiled path is 4.6x faster than the kernel chain it replaces. What this does NOT say: the scene has three sources and the standard preset has 577, so the radiation loop of a real epoch is far larger on both sides and the kernel's share of `Model.accel`'s 88 ms stays unmeasured (O-31, narrowed). **Superseded: there is no 88 ms after WP-C3, and O-31 is closed with O-29 on 2026-09-18.**

**The one failure was this session's own doing** and is recorded rather than quietly widened (rule R6, and A.4 records the same mistake once before): `tNgfMex/singleEpochSpeed` asserted an ABSOLUTE 500 µs per epoch, a number with nothing behind it, on a host that measured 901.8. Rule T5 forbids asserting host-dependent absolutes, and the test broke the rule it was written under. Decomposed: the ratio stays the assertion, the absolute becomes a smoke ceiling an order of magnitude above the measurement, so a pathological rebuild still reddens the gate and a slower machine does not.

**Not done, deliberately.** `Model` is NOT routed through the compiled path, and under P-11 there is nothing to route: `Model` is the vectorised path and `kernel.singleEpoch` is the pointwise one, both first class. What the ruling obliges instead is that they agree, which `tests/tNgfMex/singleEpochEndToEndPrecision` measures over the real arc. `bodyTorques` is not ported: it is opt-in, it is not part of a_ng, and the mirror covers it separately.

**Oracles:** the frozen Python mirror `tools/pyxcheck/ngfkernel/` (independent code, rule T9, the reason it was frozen before this level); the MATLAB kernel units themselves for the 1e-14 claim; a fault-injected copy of the C for the harness's own sensitivity.

### LP-7 — Sensors and applications ◐ (opened 2026-09-08; WP-30 delivered and GREEN, the rest outstanding)

**WP-30 `app.forceBudget` delivered 2026-09-08.** `[B, info] = forceBudget(models, epoch, r, v, q, 'Frame', 'Chunk', 'Labels')` → a table of four rows per model (drag, radiation, thermal, total) with per-axis RMS and maxima, the norm, and the share of the total. It computes no physics: it evaluates `Model.accel` in chunks and accumulates only sums of squares and running maxima, so the memory is independent of the arc length, and it stacks several models so the question a budget is actually asked — what does turning ERP or self-shadowing off cost — is one call. **Its oracle is `Model.budget`**: with one chunk covering the arc the table must be BITWISE identical, which is what pins it as a reduction rather than a second evaluation path. **Measured 2026-09-08 on the CI runner, 276/276 GREEN** (predicted 276, R2026a Update 5, 219.6 s): against `Model.budget` at one chunk **exactly 0** across all seven numeric columns; 40 epochs in chunks of 7 move the RMS by 2.9415e-16 relative and the maxima by **0**; the shares reproduce the rms_norm ratio exactly; two models stacked equal two separate calls to the last bit; explicit quaternions equal the provider's to the last bit; K = 1 gives rms == max exactly; cdScale 2 doubles the drag row and leaves the radiation row untouched, both at 0. Cost: four chunks over 20 epochs are 2.47x one chunk, i.e. the per-call overhead is real but the evaluation does vectorise over epochs (a ratio near 4 would have said it does not).

**The two dependency surfaces were finally opened** (rule: anything not exercised by a green test is a guess until the source has been read). Both findings are in Part B and both change what LP-7 can do: gjToolbox has no dynamics interface at all and `Model.rhs()`'s signature does not match the one its validator documents (B.3, O-4); lsqvce 5.3's `fit` takes three declaration kinds and declares regularisation in the data, which fixes `app.calibrate`'s shape as a batched (kind C) declaration (B.4, O-24).

**WP-29 `sensor.Electrostatic` delivered 2026-09-08.** `Electrostatic(bias, S, dr, kernel, 'Rotation', 'SampleRate', 'NoiseAsd', 'NoiseCorner')` implementing `y = bias + S R [a_ng + omegaDot x dr + omega x (omega x dr) - T dr]`, then an FIR filter. Three deliberate limits, each stated in the help rather than left to be discovered: the gravity-gradient tensor `T` is an OPTIONAL CALLER INPUT because this toolbox models no gravity (D.7 C4) and a term of order 1e-12 m/s^2 is the caller's judgement to make with a number; the CRN kernel of Wu et al. (2006) is NOT constructed, because that paper has not been read for this implementation and a guessed low-pass is worse than none — the coefficients are an argument; and the partials go through the same kernel as the observation, or they would describe a signal nobody observes. Oracle: the coupling against a twice-differentiated ANALYTIC rotation (a fixed axis with a ramping rate, so the Euler term is non-zero), and the 15 partials against central differences of `observe` itself.

**Three findings from the first WP-29 run (283 tests, 280 passed):** (i) the sign ratio is **-1.0**, so a caller wiring `referenceSystem` straight into this sensor would have had the Euler term backwards — O-32 answered and written into the help; (ii) reading the skew entries off `R' Rdot` without projecting out the symmetric part leaves the differencing error in, 6.0e-7 relative on that arc, and `referenceSystem` projects while the test's helper did not; (iii) **all three red bounds were mine, set by estimate rather than by measurement** — the third such session in a row. Two were cancellation floors (a 3e-6 bias removed from a 1e-6 reading leaves 4e-15, and 1e-15 was written) and one was a numerical-differentiation step chosen by the usual "small is better" rule when the model is EXACTLY linear and larger steps are strictly better. The fix for the metamorphic one is a bound COMPUTED from `eps` of the terms involved rather than a literal; that is the mechanism, and it should spread.

**WP-29 green 2026-09-09: 283/283, gate TRUE.** The coupling agrees with the twice-differentiated analytic rotation to 1.93e-06 against an oracle whose own truncation is 3e-7; the 15 partials against central differences to 1.96e-12 once the step exploited the exact linearity; the magnitude comparison against `referenceSystem` came out **exactly 0.0** after the symmetric projection, so the two angular velocities are bitwise negatives and the convention is settled beyond doubt. On the real GRACE-C arc the proof-mass coupling is **5.27e-10 m/s^2** — above the 1e-10 noise floor, three orders below the drag: modelled for a reason.

**WP-32 `app.densityInversion` delivered 2026-09-09.** `[T, info] = densityInversion(model, epoch, r, v, q, aMeasured, 'Frame', 'Chunk', 'MinDrag')` → a timetable of `rho_kg_m3`, the model's own `rho_model_kg_m3`, their `scale`, the projection of the measurement on the modelled drag direction, and `residualPerp_m_s2`. **The inversion is a division, not a search**, and that is physics rather than shortcut: the gas force is exactly linear in the species densities, so one model evaluation gives the whole sensitivity. Three things are deliberate and stated: it does NOT calibrate (a bias is indistinguishable from a density offset over a short arc, so laundering one into the other is refused by requiring calibrated input); it SUBTRACTS the modelled radiation and thermal recoil, so an error there goes straight into rho — which is why the perpendicular residual is reported beside every value and is asserted to move only when it should; and a NEGATIVE scale is returned rather than clipped, because clipping turns a visible failure into a plausible density floor. Oracle: the round trip (a model's own acceleration returns its own density) plus planted scales of 0.5, 1.7 and 3.25, so a unit that ignored the measurement and returned the model density could not pass.

**One red bound on the first WP-32 run (289 ran, 288 passed), and it was mine again — the fourth session running.** The claim was that inverting one epoch at a time is bitwise inverting the arc; the runner measured 1.01e-28 kg/m^3, 1e-15 relative. The inversion has no cross-epoch term, but `Model.accel` beneath it does not take the same route at K = 1 as at K = 10. **The rule extracted and written into the pitfall list: never claim bitwise across a change in BATCH SIZE** — vectorised products and norms associate by batch size. Bitwise is still claimable across a change of GROUPING that leaves the arithmetic alone, which is why the two-models-stacked and maximum-of-maxima claims in `forceBudget` hold at exactly zero. This project has now met the same mechanism three times (`pagemtimes` in `Model.evaluate`, the chunked RMS, and here) and it should not need a fourth.

**WP-19 `app.calibrate` delivered 2026-09-09, with the toolbox's own solver.** `[cal, info] = calibrate(yMeas, aModel, epoch, ...)` fitting `y = bias(t) + S a_model` with a bias basis of "constant", "linear" or "hat" (continuous piecewise linear on a uniform knot grid — local support, a partition of unity, and no end conditions; a cubic B-spline basis was NOT written, because a drifting bias is not required to have continuous curvature), an optional `S` of "none"/"diag"/"full", declared jump epochs, per-sample weights, a mask, and Huber IRLS. Two things make it small enough to own: it **separates by axis**, so it is three least squares of a few dozen unknowns rather than one of a few hundred; and it is solved by **QR on an equilibrated design matrix**, never by the normal equations — the raw columns differ by 1e7 from the units alone, so equilibration is required rather than tidy. Oracle: planted parameters recovered exactly on a noise-free arc, formal sigmas falling as 1/sqrt(K), and Huber measured against the ordinary solve on the same five gross outliers.

**Three WP-19 runs to get one test right, and the sequence is the lesson.** Rev 1 asserted the 1/sqrt(K) law on the SCALE sigma: measured 8.3 against an expected 2, because the scale's sigma is set by the spread of the modelled acceleration over the arc and not by the sample count. Rev 2 moved the claim to the bias but reached it with `'Scale', "none"`, which leaves the modelled acceleration in the residual — so sigma0 stopped measuring the noise and started measuring the signal: 0.49 against an expected 0.15. Rev 3 stopped asserting a ratio between two realisations at all and asserts the EXACT identity instead, `sigma_bias * sqrt(K) / sigma0 = 1` for a constant column, which holds per arc with no sampling spread to argue about. **The rule: separate the algebra from the statistics, and put the tight bound on the algebra.** A ratio between two noise realisations was never the right vehicle for a claim about a covariance formula.

**Second WP-19 run, 295 tests, 292 passed, and the interesting failure was a LAW I asserted that is not one.** The test claimed the formal sigma of the SCALE falls as 1/sqrt(K); the runner measured a ratio of 8.3 where 2 was expected. The reason is worth keeping: a constant bias has a constant design column, so its formal sigma is exactly sigma0/sqrt(K), but the scale's column is the modelled acceleration, and what determines its sigma is the SPREAD of that acceleration over the arc — 240 epochs of a GRACE orbit do not carry four times the spread of 60. The law is now asserted where it is a law (bias, `'Scale', "none"`, bound 0.15 relative, the residual spread being sigma0's own ~5 per cent) and the scale's ratio is RECORDED beside it rather than bounded by a wish. Two smaller ones: a stale dependency count of 4 that the lsqvce removal should have taken with it, and a singularity provocation that did not provoke — a 2 s knot interval is estimable at this sampling and merely useless, so the refusal test now uses 0.5 s, where the knot count exceeds the sample count whatever the fixture's rate, and the condition branch gets its own provocation paired with a call that must succeed.

**First WP-19 run: Code Analyzer, not a test.** Removing the lsqvce row left `localContentsVersion` with no caller, and `codeIssues` flagged DEFNU before a single test ran. Removed rather than kept "in case" — every other dependency publishes a version function, which is the better probe anyway. Worth noting as a pattern: **dropping a dependency leaves dead code behind it**, and the Code Analyzer gate is what catches that, which is an argument for it running before the suite rather than after.

**WP-29 `sensor.AtomInterferometer` delivered 2026-09-09.** `AtomInterferometer(kEff, T, cycle, axis, 'Rotation', 'SampleRate', 'PhaseBias', 'GainError', 'PhaseAsd')` giving `dphi = phaseBias + (1+g) k_eff T^2 [<n,a> - 2 <n, omega x v>]`. Three things it refuses to hide. **The phase is an average, not a sample**: the pulse sequence responds through a triangular sensitivity of half-width T, applied here as a normalised kernel whose unit sum makes a constant acceleration reproduce `k_eff a T^2` exactly — the identity the reference test asserts at 1e-15. **Dead time**: one reading per cycle and nothing interpolated between, because how to bridge the gap is a decision about the caller's problem. **The 2 pi ambiguity**: `y.value` carries the full modelled phase because a model knows its fringe, `y.info.wrapped` carries what an instrument actually returns, and `y.info.fringe` the integer between them; the reconstruction is asserted. The Coriolis term is not an optional refinement. **Corrected 2026-09-15 (LP-C F6):** this paragraph stated 6.73e7 rad, an apparent 16.7 m/s², from 1.1 mrad/s across the 7.6 km/s ORBITAL velocity — but atoms and satellite share that velocity and it cancels; the term is driven by the atoms' velocity relative to the instrument, now the input `atomVelocity`. At 1 mm/s it is 8.86 rad, an apparent 2.2e-6 m/s², seven times a 3e-7 drag. Both the phase and its ratio to the drag are recorded by the reference test.

**WP-29 `sensor.DragFree` delivered 2026-09-09 — LP-7 IS COMPLETE.** `DragFree(bandwidth_Hz, order, 'Authority', 'SampleRate', 'Bias', 'GainError', 'ThrusterAsd')` returning the residual a drag-free loop FAILS to remove: `bias + (1+g)[S{a_ng} + saturation excess]`. It is the inverse of an accelerometer — `sensor.Electrostatic` reports what a sensor reads, this reports what a controller cannot cancel — and it sits in the sensor layer because that is what a payload downstream sees. Four things stated rather than assumed: the filter is a **cascade of first-order bilinear sections, not a Butterworth**, chosen because C2 rules out the Signal Processing Toolbox and a hand-written IIR is four lines that give the right asymptotes; **saturation is modelled**, because a finite thruster authority is the failure mode a storm produces and leaving it out would make the residual look linear exactly where it stops being so; the **start-up transient is masked, not trimmed**, since it is an artefact of where the arc was cut; and the noise model follows the **complementary** sensitivity, so a wider loop rejects more drag and admits more thruster noise — the design trade is in the class rather than in a comment. Oracle: the analytic |S(f)| on the unit circle, derived from the same difference-equation coefficients by a different route, checked over three decades.

**First DragFree run: 307 tests, 306 passed.** Two notes. The predicted count was written as 308 and the arithmetic was mine: 301 minus the stub row plus seven is 307, and the number in the commit message was wrong before the run was. The one failure was a bound below its own cancellation floor for the third time in this project — recovering a 1e-9 bias by differencing two residuals of order 1e-7 cannot beat eps(1e-7) = 1.3e-23, and 1e-24 was written from an estimate. It now COMPUTES the floor as ten times `eps` of the quantity being differenced, which is the mechanism already used in `tNgfApp` and which should simply be the default for any cancellation-limited claim.

**O-34 closed 2026-09-09, and it is the last thing that stood between this repository and being public.** See Part H.

**LP-7 is closed. The four remaining stubs are blocked on data this project does not hold** (`readMAG1B` O-8, `importCERES` O-3, `igrfCoefficients` and `env.magneticField`), so they are not LP-7 debt. Next is **LP-R**.

Original prompt: WP-30 `forceBudget` (moved here from LP-6 on 2026-09-08: the application layer, a reduction over `Model.budget`) · WP-29 `sensor.Electrostatic` (bias, 3×3 scale, CoM coupling from ω, ω̇ via `angularVelocityFromRotationRate`, CRN filter; oracle: [J] §3 and the GRACE-C CoM calibration of Hu 2022), `AtomInterferometer` (Δφ = k_eff a T², dead time, 2π ambiguity, rotation input; oracle Lévèque 2022 CARIOQA numbers), `DragFree` · WP-19 `calibrate` through lsqvce (oracle: recover the published GRACE-C along-track scale ≈1.05 and bias order from real ACT1B; recover planted parameters from synthetic data to their σ) · WP-32 `densityInversion`.

### LP-R — Release ◐ (opened 2026-09-09)

**The C2 toolbox scan is delivered (part of WP-21), and it earned its place on its first run.** `auditPackage.py` check 9 refuses any call into Signal Processing, Statistics, Optimization, Curve Fitting, Symbolic, Aerospace or Mapping, across code AND help text — because a reader who copies an EXAMPLE that needs a product the toolbox refuses to need has been misled by the toolbox itself, and `rms()` had reached exactly that in `sensor.DragFree`'s example plus a test helper on the same day.

Two design points worth keeping. **A name the file assigns to is a variable, not a call**: `Model.m` and `app.forceBudget` both keep a matrix called `rms` and index it, and flagging those would have taught everyone to ignore the scan on its first run — so each file's assignment targets are collected and excluded within that file. **The list is deliberately conservative**: `fit`, `smooth`, `filter` and `median` are left out because they collide with ordinary identifiers, and a check people learn to work around is worse than no check. A line that only NAMES one of these carries `toolboxScan:allow`, so every exception is visible in the diff that introduced it. Fault-injected with a `filtfilt(` call and it turned red.

Check 8 changed with O-34: there is no vendored source left to digest, so it now refuses only a PARTIAL `mex/nrlmsise00` tree — the one state that would compile into something nobody chose.

**WP-21's documentation build is delivered (2026-09-09), and O-11 is closed.** `tools/buildDocs.py` renders every page plus `helptoc.xml`; the audit compares by content and was fault-injected. `buildFunctionList.py` is folded into it, so one generator owns `doc/html` (R1). The generator caught its own reason for existing on the first run: `Sensors.html` described a CRN filter that `sensor.Electrostatic` deliberately does not construct, and called three delivered classes future work.

**WP-22 demos delivered 2026-09-09, and the gate RUNS them.** Four scripts under `demos/`, under the names `demos/Contents.m` had planned for them: `demoForceBudget` (two fidelity tiers side by side on the real arc), `demoDensityInversion` (a planted scale recovered), `demoGraceFOCalibration` (a planted bias and scale recovered — **no ACT1B file exists**, on the host or here, O-8, so planting is the only honest option and the demo says so) and `demoCAI` (the cold-atom scenario, whose headline is the 6.7e7 rad Coriolis phase rather than the drag). `tests/tNgfDemos.m` discovers them from the folder and executes each one, asserting completion and that nothing NaN or Inf is presented as a result — **not** the numbers, which are already bounded by the tests of the units the demo calls, and asserting them twice would put a second authority in a file whose purpose is illustration. What it catches is the failure a demo actually suffers: an API that moved. The reason for running them at all is on the record: a demo nobody executes rots exactly as `doc/html/Sensors.html` did. The C2 scan now covers `demos/` too, and `setupNonGravForces` puts the folder on the path.

**WP-22's benchmark and packaging delivered 2026-09-09; LP-R has nothing left but the release itself.**

`tools/benchmark.m` times the paths a user waits on and writes them to `tools/benchmark_last.txt` (appended, gitignored, the human log) and `tools/benchmark_reference.csv` (overwritten, committed, the authority) with the host and MATLAB release beside them. It repeats no unit's `speed` test: those are RATIOS by design, because a ratio survives a change of host (T5). What no test produces is the per-epoch cost, which is host-dependent by nature and has no business gating anything — **and which is exactly the figure this project quotes most, the 88 ms of O-29, measured once by hand and repeated ever since. A number that cannot be regenerated is a number nobody can check.** It is deliberately NOT run by the gate: it measures the machine as much as the code, and a gate that fails because a laptop was busy teaches people to ignore gates.

**It was first run on 2026-09-18, nine days after it was written, and it did not work.** The native row placed its table span at `[60, 60 + dt]`, which with a 120 s margin asks for attitude from -60 s while the fixture's attitude starts at 0, so `prepare()` refused and the whole function errored for EVERY value of `'Epochs'`, its own default included. It had been written, committed, cited in `Model.m`'s help and described in this document, and never once executed. The span is now centred in the attitude the SCA1B record actually carries, measured from the record; the refusal is caught and REPORTED as a row rather than pre-empted, because re-deriving `buildTables`' outward grid snap here would put that rule in two places and compare them never. **And once it ran, it was still measuring wrong**: a single warm-up call is not warm, and the first properly warmed run moved `kernel.singleEpoch` from 2.44 ms to 0.456 and the fast preset at K=1 from 31.9 ms to 4.00 — factors of 5.4 and 7.8 of pure measurement artefact. An argument had already been written on the artefact (that the presets differ by 4 % at K=1, so K=1 measures no physics); the truth is that they differ by 4.0x there, and it was withdrawn. It warms five times now. `'Epochs'` also became a VECTOR, defaulting to [1 10 100], because a per-epoch cost quoted at one batch size is what started all of this. Its output is committed at `tools/benchmark_reference.csv` and `auditPackage.py` check 13 now compares every per-epoch figure quoted in the repository against it. **A number nobody can check was the diagnosis; an instrument nobody runs is the same disease; and an instrument run once, badly, is how the cure nearly shipped its own defect.**

`tools/packageToolbox.m` was delivered and then **removed the same day**: the owner ruled that the toolbox is published via git only, so a packaging script would be a second release path nobody runs and everybody believes (O-36). `tools/benchmark.m` lives in `tools/` rather than the package — a development instrument is not part of the scientific surface, and putting them in `+NonGravForces` would have bought a catalogue row and six meaningless tests (a metamorphic test of a stopwatch) in exchange for nothing.

**What the package deliberately leaves out:** `mex/nrlmsise00` (never ours, O-34), the built MEX binaries (per-machine), `tests/` and `tools/`, and `data/`. That last one has a consequence worth stating rather than discovering: **the demos do not run from an installed copy**, because they need the 2020-01-01 GRACE-C arc. That was O-35, and the owner answered it the same day: **data is optionally downloaded during setup, and otherwise the user supplies it.** `setupNonGravForces` now fetches the IERS and CelesTrak files when they are missing; GRACE-FO Level-1B needs an Earthdata login and is named, not fetched.

Original prompt: WP-21 documentation build (`tools/buildDocs.m` regenerates `doc/html` from help blocks and the catalogue; `scanHelp.m` enforces the template; a static scan refuses Signal/Statistics toolbox calls) · WP-22 demos, benchmarks, `.mltbx` packaging, CI. No compatibility layer (F2).

**LP-R waits for LP-C** (owner, 2026-09-15): the release is cut from a toolbox whose presets all evaluate and whose integrator entry is usable.

### LP-C — Completion ◐ (opened 2026-09-15)

**Why this level exists.** A whole-repository review on 2026-09-15 (fifteen findings F1–F15, reported through the review tool) found two physics defects in the standard preset, an integrator entry point no integrator can call, a `precise` preset that cannot evaluate one epoch, and options that are accepted and then ignored or refused. The owner ruled the same day (Part I): fix all of it; make self-shadowing physically correct; target `gjToolbox.DynamicsSource` with a native C entry; finish the precise preset; add HWM14; drop DTM2020, JB2008, CLL and Schamberg.

**Scope frozen 2026-09-15** (Tier A: MATLAB R2026a Update 4 over the MCP bridge on the host; CI tier through ONE draft PR from branch `lpc-completion`, one commit per work package, version bumped once at the end).

**The findings, by number, and where each is discharged.** F1 gas self-shadow deletes the force on edge-on and leeward plates (WP-C1) · F2 the Sun's shadow column zeroes the Earth heating of panels facing away from the Sun (WP-C1; O-28) · F3 `Model.rhs` cannot be plugged into `gjToolbox.DynamicsSource` (WP-C9; O-4, issue #35) · F4 `erp.anisotropy "adm"` accepted and ignored (WP-C2) · F5 Options admits values the Model can never evaluate, the whole `precise` preset among them (WP-C2, C6–C8) · F6 the cold-atom Coriolis term is fed the orbital velocity (WP-C5) · F7 the compiled path is kernel-only and nothing composes one epoch (WP-C9) · F8 the category audit credits units by substring (WP-C4) · F9 two-path parity asserted only in the fast preset (WP-C9) · F10 default mass is GRACE-C for any macro model (WP-C2) · F11 the look-up table scans every face per direction (WP-C3; O-29) · F12 provenance hashed on every inner call, and the hash is not the FNV-1a it claims (WP-C3) · F13 `series` and `partials` evaluate the environment more than once (WP-C3) · F14 `calibrate` full-scale design uses masked NaN axes (WP-C5) · F15 setup treats committed excerpts as the full products and never checks downloaded content (WP-C5).

**WP-C1 Self-shadowing, stated physically** (design approved by the owner 2026-09-15). Occlusion and illumination are two facts with two owners. `geom.selfShadow` returns the UNOBSTRUCTED fraction o_k(d) of plate k toward direction d and nothing else; the facing factor (n·d)⁺ belongs to the kernels, which already apply it. For a back-facing direction the table stores the mirror value o_k(d − 2(n·d)n), a continuous extension that no consumer weights and that keeps the interpolation honest at grazing incidence.
- Radiation, per source: F_km ∝ P_m A_k (n_k·d_m)⁺ o_k(d_m). The kernel is unchanged; the defect was only upstream.
- Panel temperatures: q_k = Σ_m α_k,band E_m (n_k·d_m)⁺ o_k(d_m). `geom.thermalState` takes 'Exposed' per source (N×K×M) and `Model` passes the radiation occlusion (F2; O-28 closes).
- Gas: molecules do not arrive as a parallel beam. New unit `geom.flowShadow` computes, per plate, epoch and species, three ratios of the occlusion-weighted to the unweighted directional integrals of the drifting-Maxwellian incident flux over the plate's front hemisphere — number flux (for the re-emitted momentum), normal momentum and tangential momentum. The speed moments are closed forms; the directions are a flow-aligned product quadrature. `kernel.gsiCoefficients` exposes the incident and re-emitted parts, and `kernel.dragPanelForce` scales each by its ratio. Limits that pin it: a plate with nothing in front of it gets exactly 1 (bitwise), so the edge-on shear returns; s → ∞ reduces to the parallel-beam occlusion along v_rel, which is what `constantCd` uses. Not modelled and bounded by the oracle: re-emitted molecules striking a second plate.
- Oracles: `tools/pyxcheck/tpmc_check.py`, a standard-library test-particle Monte Carlo of the incident number and momentum fluxes on a plate screened by a wall, written from a physics statement by an agent that did not open the MATLAB (T9); and the unoccluded directional integrals, which must reproduce Sentman's incident closed forms. Claims: unoccluded ratios exactly 1; directional integrals within 1e-6 of the closed forms for s ∈ [1, 15]; the TPMC within 4 standard errors plus the quadrature's own recorded error.

**WP-C2 Honest contracts.** Options admits exactly what evaluates. `drag.density` exponential | msis00 | msis21; `drag.wind` none | corotation | hwm14; `drag.gsi` constantCd | sentman | maxwell, the last needing `drag.specularFraction` (default empty, refused at `Model` construction); `drag.accommodation` fixed | pilinski2010 | walker2014; `erp.backend` none | pointSource | knocke | ceres; `erp.anisotropy` is REMOVED (F4: Lambertian is the model, stated in ACCURACY); `thermal.model` none | steadyState | transient (`foldedIntoSrp` removed; `transient` carries its temperatures across `series` chunks serially, refuses `parfor`, and is refused by the integrator source because a stateful model breaks gjToolbox's non-monotonic call rule); `thermal.internalPower_W > 0` needs the macro model's per-panel `internalShare`; `partials.method` analytic | finiteDiff (`complexStep` removed: measured impossible). `Model` 'Mass' is a scalar in kg or an `io.massModel` struct; omitted, it is the macro model's `mass_kg` (`geom.graceFO` carries the handbook value), and otherwise refused (F10). `MacroModel`'s unread per-panel `alpha` is removed (O-23 closes by removal).

**WP-C3 Cost.** The self-shadow table becomes an equiangular cube map with O(1) bilinear look-up; `geom.lutResolution` is the number of cells per cube-face edge (default 16, ≈5.6°), and the convergence sweep is re-run over it (F11, O-29). Provenance is built only when a caller asks for `info`, and hashed in exact uint64 arithmetic (F12). `Model` evaluates the environment once for accel and torque together and once for all scalar partials (F13). Claim: the standard preset per epoch at least ten times cheaper than the 88 ms of O-29, asserted as a ratio, recorded as an absolute by `tools/benchmark.m`.

**WP-C4 Instruments.** The category audit credits a test method to the unit with the LONGEST name prefix of the method name (F8), fault-injected with `densityInversion*` against `density`. Two-path parity in every preset arrives with WP-C9 (F9).

**WP-C5 Smaller defects.** `calibrate` full-scale rows require every model axis finite (F14). `setupNonGravForces` looks for the exact target file and checks downloaded content by its header, not by the status code (F15). `sensor.AtomInterferometer`'s rotation input is the ATOM velocity relative to the instrument, renamed `atomVelocity`; the demo, the tests and this record lose the 16.7 m/s² that came from feeding it the orbital velocity (F6).

**WP-C6 NRLMSIS 2.1** (`drag.density = "msis21"`). The unmodified NRL Fortran (`nrlmsis2.1.tar.gz`, SHA-256 `41e47b29f795d36a5cc252b2858aa2a384c4a7323ace3d48d3ea2f2b37a1a6a8`, fetched by `setupNonGravForces`, never redistributed: Open Source Academic Research License, non-commercial) compiled with the gfortran 14.2 of the MinGW-w64 support package, double precision (`-DDBLE`), called from a C gateway through a bind(C) Fortran shim of this toolbox's own. Oracle: `msis2.1_test_ref_dp.txt`, NRL's 200 reference cases, to the printed digits.

**WP-C7 Accommodation** (owner ruling 2026-09-15: both open models, `precise` uses the second). `pilinski2010`: α = K n_O T / (1 + K n_O T), K = 4.98e-17 m³ K⁻¹, fitted with NRLMSISE-00 n_O and T over 150–500 km (Pilinski 2011 thesis §7.3, eqs 7.1–7.2). `walker2014`: adsorbate coverage θ = K_L P_O / (1 + K_L P_O), P_O = n_O k T in Pa, K_L = 1.44e6 Pa⁻¹, the adsorbate state fully accommodated and the clean substrate at α_s = 3μ/(1+μ)², μ the mole-weighted mean mass of H, He, N, N₂, O, O₂ over the surface particle mass (Paul, Sheridan, Licata & Mehta 2022, eqs 10–13, which attribute the coverage law to Walker et al. 2014 and α_s to Mehta et al. 2022); the two states are mixed in the per-species re-emitted momentum with weight θ. The surface particle mass is required, never defaulted; `precise` states 263.3223 amu (owner ruling). SESAM (Pilinski 2013) is NOT implemented: its sticking law and shape factor are in no source in hand (O-16 narrowed). Oracles: the published formula values at stated inputs; the metamorphic limits θ → 0 and θ → 1.

**WP-C8 CERES** (owner ruling 2026-09-15). `io.importCERES` reads EBAF-TOA Ed4.2 monthly means now (obs4MIPs CF files, anonymous, SHA-256 published by ESGF) and SYN1deg-3Hour once the owner supplies one real file (O-3 stays open for that product only). `earthRadiance("ceres")` uses the same cap tiles as Knocke, looks the Earth-fixed grid up at the tile centres through the ITRS rotation, bilinear in space and linear in time; albedo = SW↑/SW↓ at TOA applied to the instantaneous insolation, LW↑ emitted as given. Oracle: the grid's own global means and a pure-Python cap integration of the same grid.

**WP-C9 Integration** (owner rulings 2026-09-15). `NonGravForces.app.DynamicsSource < gjToolbox.DynamicsSource`: state 3 (GCRS), t in TAI seconds from `context.Epoch`, attitude from the Model's provider or the nominal attitude of the state, optional estimated `cdScale`/`crScale` with analytic parameter partials. `Model.rhs` is removed (superseded; O-4 and issue #35 close). Native entry: `mex/ngfRhs.c` evaluates the whole chain for one epoch in C over tables built in `prepare()` (R, Rdot, Sun, Moon, attitude nodes, drivers, mass, CERES grid, occlusion cube map) and publishes a `gj_native_entry` (ABI 1). Vermeille's geodetic conversion and quaternion slerp are PORTED from referenceSystem into this C and verified against referenceSystem's MATLAB (owner ruling: a checked copy, R1 with a mechanical comparison, upstreaming a later option). Oracles: `Model.accel` at the same epoch and state in every preset; gjToolbox's own conformance harness for the entry.

**WP-C9 design, refined before the code (2026-09-16), from reading odeSuite at `6a3f382` (`gjToolbox.DynamicsSource`, `csrc/gj_native.h` ABI 1, `Dynamics.m`, `propagate.m`, `VariationalSource.m`):**
- `app.DynamicsSource(model, 'Epoch', tai, 'Attitude', "provider"|"nominal", 'AlongVelocity', tf, 'Parameters', ["cdScale" "crScale"], 'Native', tf, 'Gravity', src, 'Margin_s', 3600)`. StateSize 3 (GCRS), `Vectorized = true`. It adds the non-gravitational acceleration to an optional `Gravity` source's (another `gjToolbox.DynamicsSource`, prepared first; its native entry, if it publishes one, is called by ours), because the native stepper binds ONE entry and NonGravForces owns no gravity.
- The epoch is the constructor's and, when the context carries one, must equal `context.Epoch`: `propagate` re-prepares a variational run with a context that has no Epoch and a StencilSpan equal to the span, so the tables are built over the span plus `Margin_s` on either side, and an evaluation outside them returns `GJ_NATIVE_EDOMAIN`.
- **One model from every call path.** With `Native` true (the default), `prepare` builds every time-only quantity into tables — the ITRS→GCRS rotation and its rate (cubic Hermite on R with Rdot, Lagrange on Rdot), Sun and Moon (Lagrange), the attitude nodes already sign-continued exactly as `referenceSystem.interpolateQuaternions` does (slerp in C), the density and wind drivers per 3-hour UTC bin, the UTC−TAI offset, the mass, the CERES months, the self-shadow cube map — hands them to `mex/ngfRhs*.c`, which evaluates the WHOLE chain for an n×m block, and publishes the entry. `accel()` in MATLAB evaluates that same entry through our gateway. `Model.accel` stays the normative reference; parity is asserted in every preset at table nodes and bounded between them.
- Native scope: every Options value except `geom.selfShadow "raycast"` (use `lut`), `thermal.model "transient"` (stateful; refused by the source in both paths) and attitude `squad`; refused at `prepare` with the offending option named, never degraded.
- Breakpoints: the drivers change at 3-hour UTC boundaries; the source declares them aligned to the grid as "restart", and the tables switch AT the aligned node, so declaration and model agree exactly (up to half a step from the true boundary, stated).
- Gradient: `hasGradient` is true only with `Parameters`; Gr and Gv are central differences of the same entry, Gp analytic (`Model`'s parameter partials); without parameters the integrator's own fallback is the honest answer.

**WP-C10 `env.magneticField`** (IGRF-14, the last stub but `importCERES`): Schmidt semi-normalised synthesis in ITRS and body frame. Oracle: an independent standard-library Python synthesis and published check values.

**WP-C11 HWM14** (`drag.wind = "hwm14"`, opt-in in every preset; owner ruling 2026-09-15). Unmodified Fortran and data files fetched by setup, compiled like WP-C6. Oracle: the test-program output shipped with the model.

**WP-C12 Records.** This document, `CLAUDE.md`, `THIRD-PARTY-NOTICES.md`, `doc/functionlist.csv`, `Contents.m`, the generated `doc/html`, the theory §2.3 drag/lift correction, the version.

**Delivered on Tier A, 2026-09-15 and 2026-09-16** (MATLAB R2026a Update 4, `matlab -batch` over an isolated dependency root — junctions to clean `referenceSystem` and `odeSuite` clones — because the sibling `referenceSystem` was checked out on another session's branch). **The level closes at 365 passed, 0 failed, 0 incomplete, 0 category shortfalls, predicted 365** (2026-09-16, 159.2 s, `tests/last_run.txt`), with `codeIssues` clean on `+NonGravForces`, `auditPackage.py` clean and `mex/test/runCChecks.sh` clean at 14 bounds. Rounds on the way: 334 (C1–C4, C7), 352 (C8, C10), 358 (C11), 365 (C9, C12). **GREEN ON CI TOO** (2026-09-16, run `35054822337`, both jobs): `auditPackage.py` and `runCChecks` on ubuntu, then the self-hosted Windows runner at `GATE TRUE (365 passed, 0 failed, 0 incomplete, 0 audit shortfalls)` with Code Analyzer clean. The runner reached it only after the owner widened the read token to `astrobayer/odeSuite` (O-39), and that wall had been hiding two defects this host could not see — see O-39. **LP-C is ☑ on both tiers.**
- WP-C1–C4, C7: **334 passed, 0 failed, 0 incomplete, 0 category shortfalls** (predicted 333 + 1 for `accommodationReference`).
- `flowShadow` against the TPMC oracle (18 cases): worst |z| 2.58 at quadrature [96 128]; against Sentman's incident closed forms 1.6e-6 at the default [12 16]. `accommodation` against its formula oracle (12 cases, 5 columns): 0 relative (bound 1e-13); fault-injected, a 1e-12 change of the surface mass reads 9.3e-13 and a tenfold argon density reads 0, as it must.
- **F1 measured** on 120 real epochs (2020-01-01 GRACE-C, standard preset, every fifth GNV1B sample of the committed excerpt), before = `main@8b62e04`, after = this branch: the RMS drag acceleration with self-shadowing ON was 0.753 of the value with it OFF before, and is 0.9974 after. The correction moves the drag by 1.68e-9 m/s² RMS, almost all along track (per axis 1.68e-9, 1.21e-10, 7.22e-11).
- **F2 measured** on the same epochs: the nadir panel's steady-state temperature was 0 K at 11 of 120 epochs before (its Earth heating zeroed by the Sun's shadow column), and is ≥ 244.1 K after. Thermal recoil moves by 1.77e-9 m/s² RMS (3.44e-9 max), radial. Radiation pressure moves from 2.4203e-8 to 2.4237e-8 m/s² RMS. With self-shadowing OFF, before and after agree to 0 — the change is confined to the defect.
- WP-C3 cost: 1.21 ms per epoch in the standard preset, 0.30 ms in `fast`, 2.06 s to construct (host, recorded not asserted), against O-29's 88 ms. **Both figures are superseded by the measured table in `tools/benchmark_reference.csv` (0.798 ms at K=100, 2026-09-18); what matters here is that this line discharged O-29 and the row stayed open anyway.**
- WP-C6 NRLMSIS 2.1: **334 passed, 0 failed** (predicted 334: the checks extend existing methods). `setupNonGravForces` fetched NRL's tarball, the SHA-256 matched the pin, `buildMex('Target', "msis21")` compiled NRL's seven files unmodified plus the shim in 2.0 s (MinGW-w64 gfortran 14.2, `-O2 -cpp -DDBLE`, libgfortran and libquadmath linked statically: the binary imports no compiler DLL, read from its import table). NRL's 200 reference cases: worst 0.9986 of half a unit in the last printed digit (bound 1); fault-injected — O and N₂ swapped reads 2.2e12, the ap history instead of daily Ap 1.2e4, ρ × (1 + 1e-4) reads 2.84. Missing species (103 rows) are returned as exactly 0. The model's mass density against its own species and masses: 3.3e-16; against `const.speciesMass` 3.1e-4. The anomalous-oxygen flag moves ρ by exactly n(O*) m_O to 8.5e-17. Cost 19 µs per epoch, 5.1× `msis00`. At 490 km on 2020-01-01 msis21/msis00 is 0.88–0.99.
- WP-C5: **334 passed, 0 failed** (predicted 334). F14: `calibrate` with `Scale "full"` now refuses a model NaN on any axis of an epoch another axis uses (before, it reached the normal matrix and surfaced as `numericSingular`, naming the wrong cause). F6: `atomVelocity` replaces `velocity`, the old field is refused, and the recorded Coriolis size is 8.855 rad (7.3 drags) at 1 mm/s. F15, verified end to end against the live servers on 2026-09-15: IERS, CelesTrak and the GFZ Level-1B archive fetched and content-checked in 26 s, a second call touches nothing. It found a defect nobody had seen: MATLAB's `untar` takes no member list, so the Level-1B fetch had raised `MATLAB:narginchk:tooManyInputs` on every attempt since it was written and reported it as a download failure. Setup now reads the members out of the tar itself, and writes them to `data/l1b` — never into `data/samples`, whose committed excerpts carry the same names. `ciGate` passes `'FetchData', "never"`.
- WP-C8 CERES: `io.importCERES` reads the obs4MIPs EBAF 4.2.1 files (rsut, rlut, rsdt; anonymous from ESGF, SHA-256 pinned in setup); `earthRadiance("ceres")` puts them on Knocke's cap tiles on the product's 20 km reference sphere, geodetic bilinear in space and linear in time, albedo = SWup/SWdown on the instantaneous insolation, longwave as given. A two-month excerpt (2019-12, 2020-01; CC BY 4.0) is committed as netCDF-3 so a standard-library reader can check it. Oracles: `tools/pyxcheck/ceres_check.py` (agent, T9): 48 cells bitwise and 6 global means; the Earth irradiance at 8 cases — tiling error 3.5e-3 of the total at [24 48], ratio to [48 96] 4.0, Richardson residual to the oracle's 1.7e-5 quadrature. Two facts measured on the way: EBAF's SWup carries 0.3–0.6 W m⁻² of twilight flux where SWdown is 0 (lost: albedo 0 there) and ratios up to 11.8 at the polar-night edge; and MATLAB's `regexp` returned no token from inside a `(?:...)?` group even when it matched, which would have read a CF time origin of 12:30:15 as midnight.
- WP-C10 `env.magneticField`: IGRF synthesis in tesla, pole-safe (P_n^m = s^m T_n^m), time taken from `io.igrfCoefficients` (one authority for the time rule). Oracles: `tools/pyxcheck/igrf_check.py` (agent, T9; exact rational Legendre polynomials, a different route): 12 cases to 8.7e-16; the GRACE-C raw magnetometer on 2020-01-01: intensity ratio 0.9961–0.9985. At a pole north and east are undefined by the position (the oracle names them by the case longitude), so the Cartesian vector is compared there.
- WP-C11 HWM14: `env.wind("hwm14")` subtracts NRL's HWM14 quiet-time wind and, with the current 3-hour ap, the DWM07 disturbance wind, rotated from the local geodetic frame through the ITRS to the GCRS on top of the co-rotation; `Model` passes its space weather. Source: NRL's release files, and the two data files NRL's server answers 403 for from rilma/pyHWM14 at a pinned commit — all seven SHA-256-pinned in setup, and vouched for by NRL's own check driver, built from exactly these files, reproducing `Check/gfortran.txt` line for line at -O0 and -O2. Oracle: that check output, through `env.wind` (1995 dates and referenceSystem's full EOP series), rows the unit cannot pose named and left out. **Found by running:** HWM14 never returns for a negative altitude — a geodetic round trip at the surface gives -1e-9 km and hung a batch MATLAB for 30 minutes; the unit and the gateway refuse below the ellipsoid (measured standalone: -1e-6 km and -1 km time out, 742 and 5000 km answer).
- WP-C9 Integration: `app.DynamicsSource` serves gjToolbox's frozen interface and `Model.rhs` is **gone** (O-4 and issue #35 close by removal, not by an adapter round it). `mex/ngfRhs.c` (≈1100 lines) evaluates the whole chain for an n×m block in C over the tables `prepare()` builds, publishes a `gj_native_entry` (ABI 1) that `gjToolbox.DynamicsSource.checkNativeEntry` accepts, and calls an optional gravity source's entry itself, because the native stepper binds one entry and this toolbox owns no gravity. `tests/tNgfIntegration.m` (7 methods, all seven categories). **Parity against `Model.accel`, which stays normative, in four configurations** — exponential/constantCd/no shadow, msis00/sentman/Knocke/cube map/steady, msis21/walker2014/CERES/SOLAARS-CF/lut32, and hwm14/maxwell/pointSource/lunar: worst 2.504e-14 relative AT the table nodes and 8.124e-11 between them at a 60 s step (the tables' interpolation error, not the model's). Fault-injected: a 5e-10 relative change inside `ngfRhsAccel` is caught by the node bound and by the two-path precision check, and is invisible to the metamorphic and vectorisation checks, as ratios must be. One `accel` call is 0.0234 of the MATLAB path (43×); a 300 s arc through `gjToolbox.propagate` gives **the same trajectory to 0 m** on the vectorized and the native binding, and a 480 s arc with odeSuite's own native J2 composed inside our C agrees with the MATLAB path to 0 m. Recorded, not asserted: 187 evaluations in 4.16 s (MATLAB accel) against 0.26 s (C accel + C stepper), 15.8×.
- Two traps in the NRL code, both handled without editing it and both documented in `mex/ngfMsis21Shim.F90`: MSISINIT holds the parameter path in 128 characters (the session scratch path measured 146, so the full name is split across MSISINIT's path and file arguments, admitting 256); and a missing parameter file or a NaN altitude reaches a Fortran `STOP`, which ends the MATLAB process — the shim checks the file first and the gateway refuses non-finite inputs before the model sees them.

**Out of scope, recorded:** DTM2020 and JB2008 (their SET drivers stay readable in `io.spaceWeather`); CLL and Schamberg (O-17); SESAM 2013 (O-16); CERES angular distribution models; multiple reflections of re-emitted molecules.

### LP-D — Gas-surface interaction, from the sources in hand (opened 2026-09-16)

The owner supplied two papers on 2026-09-16 that had been refused to the assistant (HTTP 403 on
2026-09-07). **Read first, scoped second.** They do not unblock what O-16 and O-17 ask for, and the
first duty of this level is to say so instead of implementing something adjacent and calling the
obligation closed.

**What the sources are.** [P10] Pilinski, Argrow & Palo, *Semiempirical Model for Satellite
Energy-Accommodation Coefficients*, J. Spacecraft Rockets 47(6), 2010, DOI 10.2514/1.49330.
[W14b] Walker, Mehta & Koller, *Different Implementations of Diffuse Reflection with Incomplete
Accommodation for Drag Coefficient Modeling*, JSR 51(5) 1522-1532, 2014, DOI 10.2514/1.A32668.
Both are filed untracked in `doc/references/`.

**WP-D1 The two Pilinski fits** (owner ruling above). [P10] eq. 8 fits **K = 7.50e-17 m^3/K** to the
Bowman & Moe sphere accommodations; the shipped `pilinski2010` carries **4.98e-17** from the 2011
thesis, so the option's name and its constant have cited different sources since 2026-09-15.
`drag.accommodation` gains `pilinski2011`; `pilinski2010` moves to the value of the paper it names.
Both the MATLAB and the C enumeration change together, and `app.DynamicsSource`'s index mapping with
them. **A domain restriction goes into ACCURACY that nothing has stated so far**: [P10] is explicit
that K was fitted with NRLMSISE-00 n_O and T and that "the application of this value should only be
performed with NRLMSISE-00" - and the `precise` preset evaluates `msis21`.

**WP-D2 Published closed-form oracles.** Until now the kernel's oracles are this project's own Python
mirror and Sentman's incident closed forms - both derived here. [W14b] eqs. 11 and 12 give the
closed-form drag coefficient of a **sphere** and of a **flat plate** under diffuse reflection with
incomplete accommodation, at stated conditions (U = 7500 m/s, T_w = 300 K, T_inf = 1100 K, s ~ 7.0).
A standard-library Python oracle evaluates them and the suite compares this toolbox's per-panel
coefficients integrated over the same geometry. That is the first oracle for the drag chain that
this project did not write.

**WP-D3 The deviation from the state of the art, quantified.** [W14b] measures what the closed-form
incomplete-accommodation treatment this toolbox implements (Sentman with the Koppenwallner
re-emission ratio, its T_k,r form) costs against a physically correct CLL reflected velocity
distribution: **~1.7% in C_D for GRACE at alpha = 0 and ~1.4% at alpha = 0.9**; for a flat plate the
difference falls from ~2.6% at normal incidence to ~1.4% at 77.5 deg. It belongs in the ACCURACY
block of `kernel.gsiCoefficients` and in theory section 2.3, because it is the answer to "is the
approximation justified" and the toolbox has never carried a number for it.

**Out of scope, and what would change that.** SESAM (O-16) needs Pilinski, Argrow & Palo,
*Semiempirical Satellite Accommodation Model for Spacecraft Aerodynamic Analysis*, JSR 50(3), 2013,
DOI 10.2514/1.A32348. Closed-form CLL (O-17) needs Walker, Mehta & Koller, *Drag Coefficient Model
Using the Cercignani-Lampis-Lord Gas-Surface Interaction Model*, JSR 51(5), 2014, DOI 10.2514/1.A32677
- the companion paper in the same issue, which is also the home of the adsorbate coverage law
`walker2014` uses. **Correction to this session's own first reading:** the REFERENCES block of
`kernel.accommodation` already disambiguated by PAGE (51(5) *1544*, the CLL paper), so the citation
was right all along; only the ATTRIBUTION prose's bare "JSR 51(5)" was ambiguous between two papers
by the same three authors in one issue, and it now carries the page. What A32668 does establish is
that it is *not* the source: it mentions the Langmuir isotherm once and credits it to Pilinski.

**Delivered on Tier A, 2026-09-16** (R2026a Update 4, `matlab -batch ciGate` over an isolated
dependency root with `NGF_GJTOOLBOX` cleared): **366 passed, 0 failed, 0 incomplete, 0 category
shortfalls**, predicted 366 (365 + `gsiCoefficientsPublishedPlate`).

- **WP-D1.** `drag.accommodation` is `fixed | pilinski2010 | pilinski2011 | walker2014`; the 2010
  option now carries the 2010 paper's own K = 7.50e-17 and the thesis refit ships as `pilinski2011`.
  The MATLAB and C enumerations moved in one commit, and the native parity test gained TWO
  configurations so that both new values are evaluated in C against `Model.accel`: worst 2.504e-14
  relative at a table node over **five** configurations. The oracle carries both fits and agrees at
  0 relative. A metamorphic check pins what no per-column comparison could: alpha/(1-alpha) is
  K n_O T, so the odds ratio of the two options is K2010/K2011 = 1.5060 whatever the state
  (measured 2.9e-15). **A domain restriction is now stated that never was**: [P10] fitted K with
  NRLMSISE-00 and says the value "should only be performed with NRLMSISE-00" - the `precise` preset
  evaluates msis21, which is outside the fit.
- **WP-D2.** `tools/pyxcheck/dria_plate_check.py` evaluates [W14b] eq. 12. At alpha = 1 their
  sqrt(T_kr/T_inf)/s and this toolbox's Doornbos-form v_re/v_in are the same number, so the
  comparison is an identity and is asserted as one: **1.6e-15 relative at speed ratios 5, 7 and 9**.
  That is the first check of the drag chain against arithmetic written outside this project.
- **WP-D3, and it is a finding, not a formality.** Below full accommodation the two formulations
  separate, because ours carries (1-alpha)/2 of the incident energy into the re-emitted speed and
  the kinetic-temperature form carries (1-alpha) 2/3. **Measured on the plate: 1.7 % at alpha 0.95,
  2.4 % at 0.90, 5.6 % at alpha 0, this toolbox reading LOWER.** [W14b] separately measures the
  kinetic-temperature form against a physically correct CLL distribution at ~1.4 % for GRACE at
  alpha 0.9. Both are of order 1-2 % of C_D at a real LEO surface, inside the 5-10 % that alpha
  itself is uncertain by, so the closed form is kept - but the toolbox no longer implies the
  deviation is zero. Recorded in `gsiCoefficients`'s ACCURACY and theory sect. 2.3.
- Two test repairs, both the test being wrong rather than a bound being tight: the precision test
  hard-coded 4.98e-17 to build the state where alpha = 1/2 (it now pins each fit at its OWN unit
  argument, so it measures the formula and not the constant), and the new parity configurations had
  to clear `drag.alpha` because `Model` REFUSES a stated alpha beside a computed accommodation
  rather than ignoring it. The second is the guard doing its job.

### LP-E — What the owner's sources unblock (opened 2026-09-17)

Four items arrived on 2026-09-17: three papers that had been refused to the assistant (HTTP 403,
2026-09-07) and two data grids from a GROOPS installation. **Read first, scoped second**, and the
first duty is again to say what they do and do not close.

**The sources.** [P13] Pilinski, Argrow, Palo & Bowman, *Semi-Empirical Satellite Accommodation
Model for Spherical and Randomly Tumbling Objects*, JSR 2013, DOI 10.2514/1.A32348 -- SESAM, the
paper O-16 names. [W14a] Walker, Mehta & Koller, *Drag Coefficient Model Using the
Cercignani-Lampis-Lord Gas-Surface Interaction Model*, JSR 51(5) 1544, DOI 10.2514/1.A32677 --
the closed-form CLL paper O-17 names, and the companion to [W14b] already in hand. [H24] Hladczuk,
van den IJssel, Kodikara, Siemes & Visser, *GRACE-FO radiation pressure modelling for accurate
density and crosswind retrieval*, Adv. Space Res. 2024. **Added 2026-09-17, the three papers this
project had been citing as one:** [M19a] March, Doornbos & Visser, *High-fidelity geometry models
for improving the consistency of CHAMP, GRACE, GOCE and Swarm thermospheric density data sets*,
Adv. Space Res. **63** 213-238, DOI 10.1016/j.asr.2018.07.009 -- GEOMETRY, not gas-surface
interaction; [M19b] March, Visser, Doornbos & Siemes, *CHAMP and GOCE thermospheric wind
characterization with improved gas-surface interactions modelling*, Adv. Space Res. **64**
1225-1242, DOI 10.1016/j.asr.2019.06.023; [M21] March, van den IJssel, Siemes, Visser, Doornbos &
Pilinski, *Gas-surface interactions modelling influence on satellite aerodynamics and thermosphere
mass density*, J. Space Weather Space Clim. **11** 54, DOI 10.1051/swsc/2021035 -- the one the
kernel's ACCURACY blocks mean and mis-cite (O-42). **Supplied by the owner 2026-09-17 in answer to
the O-6 oracle question:** [V] Vielberg & Kusche, *Extended forward and inverse modeling of
radiation pressure accelerations for LEO satellites*, J. Geodesy **94** 43 (2020), DOI
10.1007/s00190-020-01368-6. All are filed untracked in
`doc/references/`. And `data/albedo/earth_ceres_{reflectivity,emissivity}_grid2.5.dat`, two GROOPS
binary `griddedData` files.

**WP-E1 `io.importGroopsGrid`** (owner ruling 2026-09-17: the format is read in `io`, not converted
outside the toolbox). The layout was DECODED FROM THE BYTES, there being no specification in hand:
`B` + length-prefixed `griddedData` + length-prefixed version, four spaces, a uint64 flag, the
ellipsoid (a = 6378137.0, b = 6356752.31414), then count-prefixed vectors of 144 longitudes and 72
latitudes in radians, 72 row radii, 10368 solid angles, and a 10368 x 12 value table. Two checks
say the reading is right rather than merely plausible: the solid angles sum to 12.566424 against
4*pi = 12.566371, and the payload fills the file exactly.

**WP-E2 `env.earthRadiance` gains the gridded CERES climatology.** **The better file is the one the owner pointed at**: `earthRadiationPressure/CERES_SYN1deg_{albedo,longwaveFlux}_climatology.dat`, 1 degree rather than 2.5. Both are `griddedDataTimeSeries`, decoded from the bytes: ellipsoid, 360 longitudes (centres +-179.5), 180 geodetic latitudes, 180 row radii, 64800 solid angles, then **2413 time NODES** (16 bytes each: a double fraction and an int64 MJD, spacings 28-31 days, MJD 15020..88434 = 1900-01 to 2101-01) bounding **2412 monthly grids** of 64800 float64. Nodes minus one is exactly the grid count, which is what a piecewise-constant monthly mean requires. Albedo lands at 0.11-0.75, **December 2019** area-weighted mean 0.3232 against the 2.5 degree file's December 0.3286 and January 0.3259. **CORRECTION to this session's first reading**: the values are stored COLUMN-MAJOR, one whole grid at a time, as Eigen serialises a matrix. The first probe assumed row-major and reported 0.391 and 0.657, numbers that survived into this record for an hour. Three things settle it: the two files agree only under column-major; the seasonal phase is right only under column-major (maximum in December and January, minimum in August, which is northern snow); and the corrected means sit much closer to Knocke's own a0 = 0.34 and e0 = 0.68 than the wrong ones did. **CORRECTION 2026-09-17, and it changes the design: the two file pairs do not carry the same
physical quantity, and this block had merged them.** Measured through `io.importGroopsGrid`:

| file pair | container | resolution | first field | second field |
|---|---|---|---|---|
| `CERES_SYN1deg_{albedo,longwaveFlux}_climatology.dat` | `griddedDataTimeSeries`, 2412 monthly grids 1900-2101 | 1 degree | albedo, 0.11-0.75, Dec 2019 area-weighted **0.3232** | longwave **FLUX, W m^-2**, 145-303, Dec 2019 area-weighted **235.66** |
| `earth_ceres_{reflectivity,emissivity}_grid2.5.dat` | `griddedData`, 12 monthly columns | 2.5 degree | reflectivity, monthly 0.3095-0.3286, annual **0.319** | emissivity, monthly 0.6883-0.7113, annual **0.698** |

So **only the 2.5 degree pair holds Knocke's own dimensionless quantities**, and it is that pair whose
annual means are 0.319 and 0.698 against Knocke's a0 = 0.34 and e0 = 0.68. The 1 degree pair's second
file is an OUTGOING LONGWAVE FLUX: 235.66 W m^-2 is a global OLR, not an emissivity, and feeding it to
the emissivity route would be wrong by a factor of ~340. This block previously attributed 0.698 to the
1 degree pair.

The seasonal phase that settles column-major survives the correction and is now measured on the pair
that actually carries it: 2.5 degree reflectivity peaks in December and January (0.3286, 0.3259) and
bottoms in August and September (0.3098, 0.3095), which is northern snow; emissivity peaks in July
(0.7113) and bottoms in December (0.6883), which is the northern summer surface. **It is NOT SYN1deg and O-3 stays open**:
there is no diurnal cycle in a monthly climatology. What it adds beyond fidelity is that GROOPS
uses these very grids, so a future comparison against it (O-19) becomes apples to apples.

**SCOPE FROZEN 2026-09-17, after the corrections above.** `env.earthRadiance` gains the backend
`"ceresClimatology"`, fed the **1 degree pair**, because it is four times the resolution AND because a
measured outgoing flux is strictly better than an emissivity multiplied by a mean insolation. It
therefore joins the EXISTING measured-exitance route: `.albedo` from the first file and
`.exitanceLw_W_m2` from the second, `.emissivity` left empty, and `env.radiationSources` needs no
change at all -- it already branches on whether `exitanceLw_W_m2` is empty.

- **The reference level is MEASURED, not declared.** The file carries one radius per latitude row, and
  those radii are the WGS84 ellipsoid radius at the row's geodetic latitude **plus 30 km**, to 0.02 m
  over all 180 rows. So GROOPS places these fluxes at a 30 km TOA level, and the backend reads that
  height out of the file rather than carrying a constant. It also confirms, rather than assumes, that
  the file's latitudes are geodetic -- the reader had said so on inference.
- **The latitude spacing is not exactly 1 degree** (centres 89.4966, 88.4899, ...), so nothing may
  compute a latitude index arithmetically; the axes are the file's own.
- **Time is PIECEWISE CONSTANT**, not interpolated: 2413 nodes bound 2412 grids and the container is
  defined that way. The existing `ceres` backend interpolates linearly between monthly nodes; this one
  must not, and a test asserts the difference rather than trusting the comment.
- **The native path REFUSES it.** `app.DynamicsSource` builds tables once per run and a 64800-cell
  monthly climatology is not table material; refusing is R4 and is tested, where a silent fallback to
  Knocke would be the worst outcome available. Recorded as an obligation rather than a silence.
- Tests write their own GROOPS fixture, as WP-E1's do, because the real grids are 1.25 GB and live in
  DATAPOOL; the gate must not depend on a file the runner may not hold.

**WP-E2 DELIVERED 2026-09-17.** `erp.backend = "ceresClimatology"` evaluates the 1 degree GROOPS
CERES monthly climatology on the cap tiles. Predicted 379 and ran 379 GREEN on Tier A.

**Insisting the reference level be MEASURED is what caught a wrong latitude convention, and that is
the finding.** The design said the backend must read its TOA height out of the file rather than
carry a constant. Comparing the file's 180 per-row radii against the WGS84 ellipsoid radius at the
stored latitude read as GEODETIC spreads by **71.6 m**; read as GEOCENTRIC it collapses to **0.168
m** about exactly 30000.1 m. Two independent checks agree: atan((1 - e^2) tan(89.5 deg)) = 89.49664,
which is the file's first centre to 1e-7 deg, and the spacing then runs 0.9933 to 1.0067 instead of
being uniform. **The grid is laid out uniformly in GEODETIC latitude and the file stores the
GEOCENTRIC latitude of those same points**, and `io.importGroopsGrid` had documented them as
geodetic since WP-E1. A consumer that converts before the lookup -- which is exactly what the
`ceres` backend correctly does for ITS product -- misplaces every tile by up to 0.19 deg, a fifth of
a cell. Writing `R = ell.a + 30e3` from the record would have shipped it invisibly.

**Three further corrections to this record's own reading, all committed before any code**: the 1
degree pair's second file is an outgoing longwave FLUX (145-303 W m^-2, December 2019 area-weighted
235.66), not an emissivity -- the 0.319 and 0.698 belong to the 2.5 degree reflectivity/emissivity
pair, and routing a flux through the emissivity branch would have been wrong by ~340x; a grid was
labelled January when the seek for [58849 58850] returns grids 1440 and 1441 and grid 1440 spans
2019-12-01 to 2020-01-01; and the longitude axis carries SINGLE-PRECISION jitter (1.37e-05 deg
between neighbouring gaps, first centre -179.499994, 3.2e-08 relative), so it is uniform in intent
and not in its bytes. It is therefore indexed with the AVERAGE spacing, so the jitter cannot
accumulate across 360 cells, and its uniformity is asserted at 1e-3 rather than at float noise.

**A contract fix in `io.importGroopsGrid`, found by using it.** After a `TimeRange` seek `.mjd` kept
the file's whole node list beside a seeked subset of `.values` -- 2413 nodes beside 2 grids -- so the
struct could not say which months it held and the invariant its own reference test asserts (nodes
are one more than grids) was false for every seek. `.mjd` is now trimmed to the nodes bounding the
grids returned, and the file's span moved to `aux.mjdFile`, where it is a fact about the FILE rather
than about these values. **R2: metadata travels with the data it describes, and the way to find out
that it does not is to consume it.**

**A shape defect that existed only at K = 1.** The latitude weight indexes the axis, `la(i0)`, and
for a VECTOR source MATLAB returns the SOURCE's orientation: `la` is a row, `i0` is MxK, so `la(i0)`
came back 1x(M*K) and subtracting an MxK latitude broadcast to (M*K)x(M*K). At K = 3 the shapes
still disagreed and the error surfaced; at K = 1 it produced a 72x72 array where 72x1 was meant --
the right size in one dimension. The batched path was bitwise unchanged by the fix (difference
exactly 0), so **the whole defect lived in the single-epoch call, which is the one an integrator
makes**. This is the mirror of the pitfall this project already records about batch size: there the
warning is not to claim bitwise ACROSS batch sizes, here it is that a shape can be wrong at exactly
one of them. The vectorization category is what caught it, on its first run.

**What it does NOT add.** It is a monthly climatology: the file is named SYN1deg but a monthly mean
has no diurnal cycle, and O-3 (the 3-hourly product) stays open and is a different thing. The native
path REFUSES it (O-44) rather than letting `find()` return empty and the enum silently become 0,
which is `none` -- a Model that computes Earth radiation would have become one that does not, on the
native path only.

**WP-E3 SESAM** (`drag.accommodation = "sesam"`), from [P13]. O-16 closes or the record says
exactly which constant is still missing.

**WP-E4 CLL** (`drag.gsi = "cll"`), from [W14a]'s modified closed-form solutions, whose stated
validity is ~0.5 % over the global parameter space. O-17 closes, and the same paper settles the
attribution of the adsorbate coverage law that `walker2014` currently credits to "JSR 51(5)" on
Paul et al.'s authority rather than this project's reading.

**WP-E5 DELIVERED 2026-09-17, and it begins by correcting this block's own reading of [H24].**
Both figures this scope quoted were misread, and the misreading pointed the work package at the
wrong claim.

**The 4.22 nm/s^2 is case (a), which has NO THERMAL MODEL AT ALL.** [H24] Table 3 is a ladder of
cross-track residual RMS by modelling choice, with the calibration parameters estimated from case
(e) and held fixed across all of them so the steps isolate the radiation-pressure modelling:

| case | thermal model | geometry | surface properties | RMS [nm/s^2] |
|---|---|---|---|---|
| (a) | **none** | panel | Level-1 Handbook | **4.22** |
| (b) | **instantaneous reradiation** | panel | Level-1 Handbook | **3.04** |
| (c) | thermal inertia | panel | Level-1 Handbook | 1.70 |
| (d) | thermal inertia | ray-tracing | Siemes et al. (2023) | 0.88 |
| (e) | thermal inertia | ray-tracing | fine-tuned | 0.55 |

This toolbox's `thermal.model` none / steadyState / transient map onto (a) / (b) / (c), and its
geometry and coefficients are panel and Handbook throughout. So **the configuration this block
called "4.22" is case (b) at 3.04**, and (d) and (e) are out of reach because they need a ray-traced
geometry this toolbox does not have.

**And the 0.148 nm/s^2 is not an oracle for anything this toolbox computes.** It is the RMS
difference between accelerations derived from MODELLED panel temperatures and from MEASURED
thermistor temperatures, for the SOLAR ARRAYS only, at a fitted heat capacity C = 8200 J/K and
conductivity k = 0.4 W/K. The paper bounds the reachable accuracy at ~0.6 nm/s^2 in the same
section, because two thermistors on one array disagree by that much and the model assumes a uniform
temperature per panel. **It measures their thermal model against their own thermistors, not thermal
recoil against truth.**

**So O-6 does not close, and the re-entry condition is now specific rather than hopeful.** [H24]
publishes no per-axis TRP acceleration. What it does publish is better than nothing and is now
gated:

- **A magnitude check for THIS satellite.** The residual left by omitting a term bounds the size of
  the term omitted, so case (a)'s 4.22 nm/s^2 is a GRACE-FO figure for the thermal recoil where
  Siemes 2023's 6.4 nm/s^2 is GRACE-A 2008. Measured here on the reference arc: **cross-track RMS
  5.263 nm/s^2**, which the two published figures BRACKET. `tNgfModel/modelReference` asserts both.
- **What thermal inertia is worth.** (b) -> (c) buys 1.34 nm/s^2 of residual; the same switch moves
  this toolbox's thermal acceleration by **0.950 nm/s^2** as a model difference. Those are different
  quantities -- a residual reduction against a model change -- so the assertion is that the option
  does something and that the something is the right order, not that the numbers match.
  `tNgfModel/modelMetamorphic` asserts it, with a floor that a transient model agreeing with
  steadyState to rounding would be a stub rather than a model.

**A stale measurement fell out of the same test.** `modelReference`'s comment said "measured median
7.9e-9 over this arc" while `tests/last_run.txt` had been recording 9.25e-10 for
|median - 6.4e-9|, i.e. a median of 5.475e-9, ever since LP-C's flowShadow correction moved the
thermal recoil by 1.77e-9. **A measurement copied into prose beside the instrument that records it
is a second authority for one fact**, and this is the fourth instance this month.

**WP-E3 and WP-E4 investigated 2026-09-17; NOTHING CHANGED, by owner ruling.** Reading both papers
answered four questions and opened a bigger one.

- **O-17 does not close.** [W14a]'s closed forms are for WHOLE BODIES -- sphere, cylinder, and a flat
  plate **at normal incidence only**; it defers arbitrary orientation to Storch [14], which is not in
  hand -- and its fitted parameters are both species- and geometry-specific (Tables 3 and 4). This
  toolbox's kernel needs per-panel C_p and C_tau at arbitrary incidence, which no equation in that
  paper provides. `drag.gsi = "cll"` stays refused.
- **Three assumptions this project had flagged as unverified are now verified from the primary
  source.** [W14a]'s nomenclature prints `K = Langmuir adsorbate constant, **Pa^-1**`, settling the
  unit `kernel.accommodation` had taken on inference; its eq. (31) IS the coverage law, so the
  attribution to JSR 51(5) *1544* was right and no longer rests on Paul et al.'s authority; and its
  Table 6 gives **1.44E+06** for diffuse reflection with incomplete accommodation, which is the
  model this toolbox pairs it with.
- **One shipped constant is contradicted by BOTH primary sources.** [W14a] eq. (33) and [P13]
  eq. (19) each give **alpha_s = 2.4 mu/(1+mu)^2**, citing Goodman; this toolbox ships **3.0** from
  Paul et al. (2022) eq. 13, a secondary statement of Mehta et al. (2022). [P13] also uses
  **m_s = 65 amu**, not the 263.3223 of the 2026-09-15 ruling.
- **SESAM's sticking law is in hand after all** ([P13] eqs 17-22: K_L becomes a function of incident
  kinetic energy through an initial sticking coefficient, with a binding energy near 5.7 eV), so
  O-16 is implementable. But it addresses the HIGH-ENERGY eccentric-orbit regime, where adsorption
  falls off; it would not raise the accommodation at GRACE-FO's circular 7.6 km/s.

**THE FINDING THAT MATTERS, and it is structural rather than a constant.** Measured on the reference
day (2020-01-01, 514-523 km, deep solar minimum, P_O = 1.0-2.0e-8 Pa, coverage theta = 0.014-0.028 --
three orders of magnitude below the ~1e-5 Pa where [W14a] fig. 17 says its isotherm fits):

| parameter set | alpha_s | effective alpha |
|---|---|---|
| [P13] SESAM's own: 2.4, m_s 65 | 0.272-0.343 | 0.282-0.361 |
| [W14a] constant 2.4, Mehta mass 263.3223 | 0.083-0.112 | 0.095-0.137 |
| Paul constant 3.0, Pilinski mass 65 | 0.340-0.428 | 0.349-0.444 |
| **Paul 3.0 + Mehta 263.3223 -- what ships** | 0.103-0.140 | **0.116-0.164** |
| fitted from orbit at 270-500 km (March et al. 2019, as this toolbox cites it) | | **0.82-0.93** |

**No combination of the published constants reaches the fitted values.** The surface mass dominates
-- 65 against 263.3223 amu moves alpha by a factor 2.6, where 2.4 against 3.0 moves it by 1.25 --
and even the most generous pairing reaches 0.44. At this altitude theta is so small that all three
models report essentially the CLEAN-SURFACE accommodation, and that is 0.3-0.4 at best. So the gap
is not a constant anyone has mis-copied: either the fitted 0.82-0.93 do not transfer to 520 km in
deep solar minimum, or the adsorption models mis-extrapolate far below their fitted pressure. **That
question cannot be settled from the papers in hand**: March et al. (2019, JSWSC, open access) is a
citation in this toolbox's own ACCURACY block and has not been read. Recorded as O-42.

**O-42 INVESTIGATED 2026-09-17 with the three March papers in hand. STILL NOTHING CHANGED, by owner
ruling. The implementation is NOT defective, and the question that remains is a design ruling, not a
bug hunt.** The decisive measurement was to stop comparing our 520 km numbers against a figure fitted
somewhere else, and instead evaluate OUR models at the conditions the literature actually fitted.
[M21] reports that at CHAMP altitudes under low solar activity SESAM's accommodation "oscillates
between the two fixed values of 0.85 and 0.93". Evaluated at 2008-01-01, equator, noon, F10.7 = 70,
Ap = 4 (`msis00`, this toolbox's own units, no fitting of any kind):

| altitude | P_O [Pa] | theta | `walker2014` | `pilinski2010` | `pilinski2011` |
|---|---|---|---|---|---|
| 300 km | 3.33e-06 | 0.827 | 0.857 | **0.948** | **0.923** |
| 350 km | 1.15e-06 | 0.624 | 0.686 | **0.862** | 0.806 |
| 400 km | 4.06e-07 | 0.369 | 0.468 | 0.688 | 0.595 |
| 450 km | 1.46e-07 | 0.173 | 0.295 | 0.442 | 0.344 |
| 500 km | 5.30e-08 | 0.071 | 0.194 | 0.223 | 0.160 |
| 520 km | 3.55e-08 | 0.049 | 0.167 | 0.162 | 0.114 |

**At CHAMP's altitude and CHAMP's solar conditions, `pilinski2010` gives 0.86-0.95 -- straddling the
0.85-0.93 [M21] reports for SESAM there.** The same code, unchanged, gives 0.11-0.17 at GRACE-FO's
520 km in deep solar minimum. So the factor of five is not a mis-copied constant and not a porting
error: it is what the adsorption family PREDICTS, and the family reproduces the published range where
the published range was measured. **The earlier framing of O-42 compared a 520 km deep-minimum
prediction against a 300-400 km fitted constant and called the difference a defect. It is not.**

Repeating the sweep at F10.7 = 180, Ap = 9 gives 0.63-0.83 at the SAME 520 km. So these models make a
strong, falsifiable claim this project had not noticed it was shipping: **GRACE-FO's energy
accommodation swings by a factor of five over the solar cycle**, from ~0.15 at minimum to ~0.8 at
maximum, because the atomic-oxygen partial pressure that drives the Langmuir coverage does. A constant
alpha cannot represent that, and 2020-01-01 -- the reference day this whole toolbox is validated on --
sits at the extreme low end of it.

**What the choice is worth, in the quantity the toolbox exists to produce.** Same `standard` preset,
same day, same geometry, 30 epochs, only `drag.accommodation` differing:

| accommodation | alpha | drag RMS [m/s^2] | vs fixed 0.85 |
|---|---|---|---|
| fixed 0.85 -- `standard` ships this, and [M21] recommends it | 0.850 | 6.8826e-09 | 1.0000 |
| fixed 0.93 -- Doornbos 2011, superseded by [M21] | 0.930 | 6.6158e-09 | 0.9612 |
| `walker2014` -- **`precise` ships this** | 0.116-0.164 | 8.0886e-09 | **1.1752** |
| `pilinski2010` | 0.051-0.098 | 8.1882e-09 | 1.1897 |
| `pilinski2011` | 0.034-0.067 | 8.2164e-09 | 1.1938 |

**The `precise` preset therefore predicts 17.5 % more drag than `standard` on the reference day, and
the entire difference is the accommodation coefficient.** Switching 0.85 to 0.93 moves it by 3.9 %.
For accelerometer calibration this is not a refinement, it is the leading systematic: a 17.5 % drag
error is absorbed by the estimated scale factor and propagates straight into any density inverted from
it. Two presets of the same toolbox disagreeing by that much about the same day is a defensible
scientific position ONLY if the toolbox says so out loud, which it does not.

**So the physics is settled and the ruling is not.** The open question is no longer "is our
accommodation wrong" but "which of two defensible things should `precise` mean": the state of the art
for THIS application ([M21], fitting a constant per mission, because the composition data a computed
model needs is not measured on these satellites), or the state of the art in gas-surface physics
(computing alpha from the composition a density model supplies, which is what `precise` does, and
which reproduces the literature where the literature was fitted). Both are honest. Recorded for the
owner as O-42.

**WP-E7: OWNER RULING 2026-09-17 -- `precise` takes the FITTED alpha = 0.85, and the adsorption models
become opt-ins.** The owner chose the state of the art for THIS application over the state of the art
in gas-surface physics, which is what [M21] itself does and for the reason it gives: the in-situ
composition a computed model needs "is not measured by any of the above-mentioned satellites". All
three presets now agree on accommodation, and the 17.5 % is no longer a difference between two presets
of the same toolbox -- it is a difference between a preset and a named option, which is where a
disagreement of that size belongs.

What the ruling moved, and one thing it nearly broke:

- `Options` `precise`: `accommodation` "walker2014" -> "fixed", `alpha` [] -> 0.85, `surfaceMass_amu`
  263.3223 -> [] (a value left beside a model that does not read it is the contradiction `Model`
  already refuses in the other direction).
- **`tNgfIntegration` gained a SIXTH parity configuration, and this is the part worth reading.**
  walker2014 is accommodation enum value 3 in the C evaluator, and it was reached by exactly one
  configuration: the `precise` preset. Changing the preset would have dropped it from the parity claim
  silently -- both paths would still run and would quietly evaluate different models, which is the
  identical blindness the `hwm14` case was added to prevent in LP-D. **A preset is not a test, and a
  branch covered only by a default is covered by nothing once the default moves.**
- The citation is corrected in `kernel.accommodation` and `kernel.gsiCoefficients`: "J. Space Weather
  Space Clim. 9, A36" resolved to no paper and is now [M21], JSWSC 11, 54, DOI 10.1051/swsc/2021035.
- `Model`'s refusal message for a missing `surfaceMass_amu` pointed at "the precise preset sets
  263.3223", which stopped being true in the same commit; it now names the disagreeing sources instead
  of a preset. **An error message that cites a default is a copy of that default.**
- The theory document's Energy accommodation section is rewritten around the measurement, and its
  summary table said "SESAM in `precise`", which was never true -- SESAM has never been implemented.

**Two more stale facts fell out of the same file, both four levels old, both because nothing compared
them with anything.** The theory document said "51 units" (LP-5's count; there are 56) and listed
`rhs` among `Model`'s methods (removed in LP-C). The count is now `\ngfunitcount`, generated by
`tools/buildDocs.py` from `doc/functionlist.csv` and refused when stale, like every other derivable
fact on that page. This is the fourth instance this month of the same defect and the second in the
same document.

**The rehearsal cost three runs, and none of the failures were in this repository.** Recorded
because the misreading was reasonable each time and wrong each time.

- Run 1 used `matlab -batch ciGate` with NO arguments. `DependencyRoot` defaults to the parent
  folder, so it installed from the SIBLING clones -- including the `referenceSystem` checkout that
  another session keeps on `fix/reviewFindings`, which this project is under standing instruction
  never to use. Nothing tracked was touched (that clone is still clean and still on its branch), but
  the run did `savepath`. Recorded as O-43.
- Run 2 used the isolated root correctly and reported **93 failures across whole classes**. The
  temptation was to read that as a regression from the preset change. It was not: run 1's MATLAB was
  STILL ALIVE, because the launch used both the harness's background flag and a shell `&`, so the
  "finished" notification described the shell exiting, not MATLAB. Two MATLABs then fought over one
  mapped DLL. **A whole test class failing in setup is a configuration symptom, and the diagnosis
  order -- configuration, then criterion, then code -- exists for exactly this.**
- Run 3, with nothing else running, was **374 passed, 0 failed, 0 incomplete, 0 shortfalls**, which
  was the prediction.

**A check was written for this and deliberately NOT shipped.** The obvious fix looked like promoting
`referenceSystem:setupReferenceSystem:mexBuildFailed` to an error inside `ciGate`, so the gate would
refuse instead of measuring 374 things against a binary that might be stale. It was written, and then
the green run disproved it: **that warning fires on this host every time, including in the run that
went 374/374.** referenceSystem reads its leap-second table through the gateway before it decides to
relink, so once the DLL is mapped the relink cannot succeed -- and when the binary already matches
its sources (both clones carry source hash 109793179) the failed relink changes nothing. The warning
therefore cannot distinguish "locked but current" from "stale", and promoting it would have turned a
healthy tree red. **A check that cannot fail selectively is not a check, and this one was caught only
because the green run happened to come first.** The discrimination has to live where the source hash
is known, which is upstream; recorded against O-43.

**And a methodological finding about this session's own instruments.** The first fault injection of
the new macro reported the audit CLEAN, and the obvious reading was that the staleness check does not
cover `generated.tex`. It does. The injection was a `sed` whose backslash never matched, so the file
was never modified and the audit was correctly silent about a healthy tree. **A fault injection that
does not inject proves nothing, and it fails in the direction that looks like a finding.** Assert that
the corruption landed before believing the instrument's silence -- injected with Python and verified
by reading the file back, the check fires on the first try.

**And an instrument that already knew.** `kernel.accommodation` computes `info.insideFit` -- whether
alpha falls in the 0.85-1.00 its constants were fitted over -- and NOTHING reads it; one test asserts
its size. The toolbox has been able to say "this is an extrapolation" since 2026-09-07 and never has.

**WP-E3 SCOPE, FROZEN 2026-09-17 before any code.** O-16 closes or this block says why not.

**The fact that decides the whole work package, and that would have been silently wrong if assumed:**
SESAM's P_O is NOT the thermodynamic partial pressure this toolbox already computes. [P13] eq. (15)
and Appendix B define it as a RAM IMPACT pressure, `P_O = 1/2 rho_O V_r^2 C_D,in,O(s)`, where
`C_D,in,O` is the free-molecular sphere drag coefficient with the re-emission term set to zero. At
520 km that is 4.97e-06 Pa against 3.55e-08 Pa for `n_O kB T` -- a factor of 140. Feeding SESAM the
pressure `walker2014` uses would have under-predicted its coverage by that factor and produced a
plausible, wrong, silently low accommodation.

**And the converse was checked before concluding anything about the models already shipping.**
[W14a] fits its isotherm to the THERMODYNAMIC pressure: it states its models fit theta well "for P_O
above ~1e-5 Pa (corresponding to altitudes below ~300 km)" and that "P_O is dependent on T_inf and
n_O". This toolbox measures 3.33e-06 Pa at 300 km from `n_O kB T`, which matches that statement, and
5e-04 Pa from the ram form, which does not. So `walker2014` and `pilinski2010/2011` are fed the
pressure their constants were fitted against, and the two parameterisations differ by roughly the
factor their constants differ by (K_L,o 5e6 torr^-1 = 3.75e4 Pa^-1 against Walker's 1.44e6 Pa^-1).
**Neither model is defective; they are different parameterisations and each needs its own pressure.**

The equations, all from [P13]: (17) theta = K_L P_O / (1 + K_L P_O); (18) alpha = (1 - theta) alpha_s
+ theta; (19) alpha_s = **2.4** mu/(1+mu)^2 -- note SESAM's own coefficient, not the 3.0 `walker2014`
ships from Paul et al.; (20) mu = mbar/m_s; (22) K_L(E) = s_o(E, Eb, Ta) K_L,o + K_L,f; (26b) the
closed-form initial sticking coefficient. Table C1: Eb = 5.7 eV, K_L,o = 5e6 torr^-1, K_L,f = 3e4
torr^-1, m_s = 65 amu.

**Two numerical rewrites are mandatory and are the reason this is not a transcription.** (26b) as
printed evaluates `exp(Eb/kT)` and `exp(2 sqrt(Eb Er)/kT)`, which at GRACE-FO are 1e36 and 1e65,
multiplied by an `exp(-(Eb+Er)/kT)` of 1e-66. Combining the exponents analytically first gives
`exp(-Er/kT) - exp(-(sqrt(Eb) - sqrt(Er))^2/kT)`, both in [0, 1], with no large intermediate at all.
Second, the printed `erf((sqrt(Eb)-sqrt(Er))/sqrt(kT)) + erf(sqrt(Er/kT))` cancels to nothing when
the incident energy exceeds the binding energy; the identity `erf(a) + erf(b) = erfc(-a) - erfc(b)`
removes the +-1 entirely. Both are the pitfalls this project already records for `1 + erf(x)`.

**One transcription ambiguity, declared.** The OCR of eq. (15) reads `(4s^4 + 4s^2 + 1)/(2s^4)`; the
classical Schaaf-Chambre free-molecular sphere coefficient is `(4s^4 + 4s^2 - 1)/(2s^4)`, and the
paper's own derivation (Appendix B: the sphere coefficient with the re-emission term removed) is that
classical result. The minus sign is implemented. At LEO s ~ 8 the two differ by 2/(2 s^4) = 2.4e-04
relative, which is recorded here so the choice is visible rather than assumed.

**Interface consequence.** SESAM needs the FLOW, not just the gas state, so `kernel.accommodation`
gains `'RelativeSpeed_ms'` -- required by `sesam`, refused by every other model, exactly as `'Alpha'`
is `fixed`-only and `'SurfaceMass_amu'` is `walker2014`-only. The C evaluator gains `NGF_ACC_SESAM`
and `tNgfIntegration` a seventh configuration, because **a branch covered only by a default is
covered by nothing**, and this one has no default at all.

**Predicted at GRACE-FO (520 km, deep solar minimum) alpha ~ 0.46, and at CHAMP (300 km, low
activity) alpha ~ 0.96**, computed by hand from the equations above before any code existed. If the
implementation disagrees with those, the implementation is wrong, not the hand calculation.

**WP-E3 DELIVERED 2026-09-17. O-16 CLOSES.** `drag.accommodation = "sesam"` evaluates on both
paths, MATLAB and C, predicted 374 and ran 374 GREEN on Tier A.

**Every number the scope block predicted by hand, before the code existed, was met.** The frozen
prediction was alpha ~ 0.46 at 520 km in deep solar minimum and ~ 0.96 at CHAMP's 300 km under low
activity. Measured: **0.408 and 0.962**. The pieces reconcile exactly where the hand calculation was
careful -- ram pressure 4.968e-06 Pa against 4.97e-06 predicted, sticking coefficient 0.8383 against
0.836, coverage 0.1359 against 0.136 -- and differ only where it was not: the hand value used a mean
incident mass of 16 amu by eye, where the six-species mean at 520 km in deep minimum is 11.9 because
hydrogen and helium dominate the COUNT there. **The prediction was wrong about the one quantity it
guessed instead of computing, and right about all four it computed**, which is the outcome that
tells you the implementation is sound rather than merely agreeable.

**And an independent confirmation nobody arranged.** At 350 km under low activity this
implementation returns **0.899**, and [M21] reports SESAM oscillating between 0.85 and 0.93 at CHAMP
altitudes under exactly those conditions. That figure was read weeks before SESAM was implementable
and is not a parameter of anything here.

**The oracle integrates what the toolbox evaluates in closed form.** `accommodation_check.py`
computes the sticking coefficient by composite Simpson over eq. (26a)'s DEFINING INTEGRAL, not by
eq. (26b). So the reference test is not two transcriptions of one scanned formula agreeing with each
other; it is evidence that the closed form IS the integral it claims to be -- the only claim worth
making about a transcription from a scan. They agree to 1.221e-14, far inside the 1e-9 the test
allows for the quadrature's own O(h^4) error. Independently, eq. (26b) was DERIVED from eq. (26a)
by hand (substitute u = sqrt(E); the integrand becomes 2u exp(-(sqrt(Er) - u)^2/kT), elementary), and
the derivation reproduces the implemented stabilised form term for term, including the sign of every
erf. **The scan's ambiguity at eq. (15) is therefore the only transcription risk left**, it is
declared in the help and the record, and it is worth 2.4e-04.

**The case set sweeps 3 to 11.8 km/s on purpose, and a test asserts that it does.** SESAM's sticking
coefficient is ~1 below 5 km/s and ~0 above 11; a twelve-case oracle that sat at 7.6 km/s would have
checked a near-constant and would have passed with the model's entire energy dependence deleted. The
reference test refuses a case set whose sticking never falls below 1e-3 or never rises above 0.999.
**An oracle that does not traverse the behaviour is an oracle for a different model.**

**The C moved in the same commit as the MATLAB, and the seventh parity configuration earned itself
within the hour.** SESAM is accommodation enum value 4 and the only model that needs the FLOW, so
`ngfAccommodation` had to take the relative speed it already held at the call site. The first gate
came back RED at **6.5e-02 against a 1e-13 bound** -- not rounding, a different model. The cause:
**`mex/ngfRhsGateway.c` copies the setup BY FIELD NAME**, and `eV` and `torr` had been added to the C
struct and to `app.DynamicsSource`'s packing but not to the gateway that moves one into the other.
They arrived as uninitialised memory, SESAM's Langmuir constant became garbage, and the native path
evaluated something that was not SESAM. **A preset would not have caught this and neither would any
test that ran only the MATLAB path**; the configuration added one commit earlier for precisely this
reason is what turned it up.

**That fact lives in three places, so `auditPackage.py` check 10 now compares them**: every scalar
`mex/ngfRhs.h` declares in the constant block must be named in a `scalar(st, "...")` call in the
gateway. Fault-injected by deleting the `torr` copy: it names the field and says the native path
will silently evaluate something else. It turns a three-minute red gate into a one-second refusal,
which is the difference between a check and a post-mortem.

**And a stale refusal test, the same shape as the one in `tNgfCore`.** `tNgfKernel` asserted that
"sesam" is rejected by `mustBeMember`, which was true until this commit. Both were repointed at a
name that is still inadmissible. **A test that asserts a capability is absent has to be found when
the capability arrives, and grep for the option name is the only thing that finds it.**

**Measured bounds, all with room to spare.** Native vs MATLAB over seven configurations 2.504e-14
(1e-13); the oracle's eight columns 3.799e-16 (1e-13); the closed form (26b) against quadrature of
(26a) **1.221e-14** where the bound allowed 1e-9 for the quadrature's own O(h^4) error -- so the
agreement is far better than the test demanded, which is what one expects when two routes compute
the same mathematical object rather than the same code twice. The sticking coefficient hits its two
analytic limits at 1.1e-16 and exactly 0, and K_L collapses onto K_L,o + K_L,f and onto K_L,f there.

**O-42 CLOSED 2026-09-17 by owner ruling: WARN.** `kernel.accommodation` had computed
`info.insideFit` since 2026-09-07 and no caller had ever read it, so a user who selected an
adsorption model at GRACE-FO altitude received a coefficient five times below every published fit
with nothing to say so. The owner chose warn over refuse, which is the R5 reading: extrapolating a
fitted model KNOWINGLY is legitimate, and refusing would forbid it.

- `info.insideFit` is now reported by EVERY model, not only the two that had it. `fixed` and
  `sesam` are true everywhere, and for opposite reasons worth keeping straight: a stated alpha
  extrapolates nothing, while SESAM's source states no domain this project can transcribe -- [P13]
  fits over a set of objects and orbits rather than over a range of a variable. Saying so beats
  inventing a range or leaving the field quietly absent.
- The domains are each model's own and are a DIFFERENT KIND of quantity: alpha in 0.85-1.00 for
  `pilinski2010/2011`, but a PRESSURE for `walker2014` (P_O >= 1e-5 Pa, where [W14a] fig. 17 says
  its isotherm fits), so that model's flag describes the environment it was handed rather than the
  number it returned.
- **This is the second unit in the toolbox permitted to warn**, after `dependencyReport`, and a
  deliberate exception to "no printing from library code". A warning is not printing: it carries an
  identifier and the caller can silence it. It fires ONCE PER CALL, naming how many of K epochs.

**And suppressing it correctly took two goes, which is the part worth reading.** Once the warning
existed, ~25 copies of it flooded the suite log, because most of these tests evaluate adsorption
models at GRACE-FO densities on purpose. **A run log full of an EXPECTED warning is exactly where an
unexpected one hides**, so it is off by default in `NgfTestCase`'s TestMethodSetup and `loud` turns
it on where it is the subject. The first `loud` restored the state with `addTeardown`, which fires
at the END of the test method -- so the warning stayed on for the rest of that method and three
later calls went on printing it. 25 became 3, not 0: the flood had been MOVED rather than removed.
`onCleanup` inside the helper restores it when the helper returns, and the log is clean. **A
restore that happens later than the thing it protects is not a restore**, and counting the
remaining warnings is what showed the difference between 3 and 0.

**[V] READ 2026-09-17. It does NOT close O-6, and it contradicts [H24] on a default this toolbox
ships.**

**O-6 stays open, for the reason predicted before reading**: [V] is **GRACE-A, 2010**, not GRACE-FO,
and it never publishes a per-axis thermal acceleration. Thermal appears only as a DIFFERENTIAL
effect, and the one figure given (mean differences up to 7e-9 m/s^2) lumps instantaneous thermal
reradiation together with anisotropic reflection. GROOPS `SimulateAccelerometer` remains the only
candidate that would give a per-axis thermal acceleration for the right satellite from independent
code, and it is the standing ask.

**THE FINDING THAT MATTERS IS NOT THE ORACLE.** [V] and [H24] disagree, in print, about whether
modelling thermal reradiation instantaneously helps at all -- which is exactly what `standard` and
`precise` do (`thermal.model = "steadyState"`).

| source | satellite | on instantaneous reradiation |
|---|---|---|
| [H24] Table 3 | GRACE-FO | cross-track residual 4.22 -> **3.04** nm/s^2, and 1.70 with inertia: it HELPS |
| [V] sect. 5 and conclusions | GRACE-A, one year from 2010 | fits "more than **39% worse**" in cross-track; "did not lead to improvements": it HURTS |

`kernel.thermalRecoil`'s ACCURACY block had cited [V] since LP-5 -- as support for TRP being the
second-largest radiation term -- **and had never mentioned the half of it that disagrees with this
toolbox's own default**. That is the one-sided citation this session has been finding everywhere,
and it is now recorded in the unit.

**The presets do not move, and the reasons are about evidence rather than preference**: [H24] is the
right satellite, it is four years later, and its whole Section 3 repairs the cross-track calibration
that [V] itself names as a possible cause of its own result. [V] also names neglected horizontal
winds in the drag model, which this toolbox models through `env.wind` and [V] did not. A user who
holds the other position sets `thermal.model` to "none". **It is an open question in the field, and
the ACCURACY block now says so instead of quoting the agreeing side.**

**And [V] hands O-3 a justification it did not have.** Its extended Earth-radiation model is built on
HOURLY CERES SYN1deg -- exactly the product O-3 wants -- and its Table 2 grades the radiation source
in the order this project has been implementing: Knocke 1988, then EBAF albedo/emission, then EBAF
outgoing flux, then **SYN1deg outgoing flux**. It concludes that "directly using the Earth's
outgoing fluxes improves the radiation pressure force model", with RMS fits down 5 % and 27 % for
its two calibrations. Its Table 2 also grades the DISCRETISATION, 19 surface elements against a
1x1 degree grid, where this toolbox tiles a cap -- a third option neither paper tested.

**WP-E6 Records**, including the licensing question the grids raise: they are CERES-derived and
CERES data is open, but these particular files come from a GROOPS distribution and this project
does not vendor what it cannot redistribute (the O-34 lesson). Recorded as O-41 for the owner.

---

### LP-D2 — The documentation audit and the user guide (2026-09-18)

**THE AUDIT ANSWERED TWO OF THE OWNER'S THREE QUESTIONS WITH "NO".** `auditPackage` check 2 had
verified that the INPUT/NAME-VALUE/OUTPUT sections EXIST; nothing had ever compared them with what
the functions take. Three auditors now do (`tools/audit_args.py`, `audit_class.py`,
`audit_options.py`, run by check 12). What they found on their first honest run:

- **Seven name-value arguments of `Model`** -- the toolbox's single point of access -- that a caller
  can pass, that the code validates, and that the class help named nowhere: `B_body` and
  `Dipole_Am2` on both `torque` and `evaluate`, `T0_K` on `evaluate`, `Quantities` on `series`,
  `Method` on `partials`. The first two are the magnetic-torque inputs and are needed TOGETHER or
  not at all; one alone leaves the term out silently, which is now written down.
- **Four public methods documented nowhere**: `Options.toStruct`, `DynamicsSource.diagnostics`,
  `DynamicsSource.epochAt`, `MacroModel.silhouetteDefect`.
- **Three `MacroModel` name-value arguments** whose PROPERTIES were documented while the way to SET
  them was not -- `CentreOfMass`, `ProofMassOffset`, `Reference`. The property list said what the
  field held; nothing said how to fill it, and the names differ in case from the properties.
- **Four `qc.screen` options** delegated to `qc.residuals` on one line, carrying no size, type or
  units.

Every admissible OPTION VALUE was already described, in both `Options` and its consuming unit; that
third question answers yes.

**THE AUDITORS WERE WRONG BEFORE THEY WERE RIGHT, AND THAT IS THE PART TO READ.** The first run
reported nine findings of which eight were the parser's fault: it assumed two spaces after a name
where many rows use one, it treated the token `Name` as a header when `geom.boxWing` has a real
option called `Name`, it read prose inside a table ("All arrays must be compatible...") as a row,
and a later fix let a row beginning with a capital (`F`, `T`, `GM`) end the section it was inside,
which took the count to 23. **A checker that reports what is already there teaches people to ignore
it**, which is this project's own rule about check 9. Each false positive is now a comment in the
parser, and the auditors run clean.

**THE USER GUIDE.** `doc/src/NonGravForces_UserGuide.tex`, ten chapters: introduction, installation,
theory and implementation, using the toolbox, a step-by-step GRACE-type chain, the integrator
interface, worked examples, limitations, the complete API and the references. CI builds it and
uploads it as the `userGuide` artefact, three passes because the API appendix moves the page
numbers.

- **The theory is not duplicated.** `theory_body.tex` and `references.tex` were SPLIT out of the
  theory document and are `\input` by both, so the two PDFs render the same words and cannot drift.
  The split was verified to be a move: the concatenation equals the original, line for line.
- **The API appendix is GENERATED** by `tools/buildDocs.py` from `doc/functionlist.csv` and the help
  blocks -- 56 of 56 units, 4130 lines, 148 tables, 444 argument rows -- and `auditPackage` refuses
  a stale one by content, exactly as it does `generated.tex`. A hand-written API chapter is stale
  the day it ships, which is the defect this level spent a week finding in other forms.
- **Every printed output was produced by running the code**, not composed. The four demos and the
  walkthrough were executed and their diaries pasted in.

**WRITING THE WALKTHROUGH FOUND TWO DEFECTS IN ITS OWN DRAFT CODE**, which is the argument for
running a guide's examples rather than writing them: `info.Tgas_K` does not exist, and
`app.calibrate` takes `'Bias'` and `'Scale'`, not `'Parameters'`. Shipping those untested in a
beginner's guide is exactly the rot this project records for documentation nobody executes.

**AND IT TURNED UP TWO THINGS WORTH TEACHING RATHER THAN HIDING.** The GRACE-FO panel model's
silhouette defect is **0.14 m^2, not zero** -- it is an open box-wing, not a closed body, and the
first draft's label "(0 for a closed body)" would have read as a failure. And calibrating on the
200-second excerpt returns scale factors of **0.16, 0.70, 0.85** with condition numbers 57, 46 and
7: badly determined, because a scale factor is set by the SPREAD of the modelled signal over the
arc and 200 seconds of a 90-minute orbit carries none. The guide prints that beside the synthetic
demo which recovers a planted scale to 4e-16, because the contrast is the lesson: the estimator is
right and the arc is too short.

**THE REFERENCES WERE INCOMPLETE IN A WAY THAT MATTERED.** Four sources the toolbox cites IN ITS OWN
HELP were absent from the list: March et al.\ (2021), whose fitted alpha every preset now uses;
Walker, Mehta & Koller (2014), whose coverage law `walker2014` implements and whose whole-body
closed forms are the reason `cll` is refused; Paul et al.\ (2022); and Hladczuk et al.\ (2024),
which supplies the 4.22 and 3.04 nm/s^2 of the limitations chapter. Six further entries lacked
volume, pages or DOI. Every correction comes from the paper in hand, read in these sessions, not
from memory of it.

**What this host could NOT do**: compile. There is no LaTeX engine here (O-38), so the documents are
checked structurally -- environments balanced, `\input` targets present, no unsubstituted
placeholder, no unescaped underscore outside math -- and CI is the compiler. That checker also had
to be corrected: stripping only `$...$` reported 24 underscores, every one inside an `equation` or
`align` block of a document CI already builds.

### The fourth time measuring overturned something this record asserted (2026-09-18)

O-43(b) was written on 2026-09-17 as an upstream defect: `setupReferenceSystem` warns
`mexBuildFailed` whether a binary is stale or merely locked, a consumer cannot distinguish, file it
upstream. The owner then asked for it to be filed. **Measuring first found the premise wrong.**

- That setup does not guess. It records a digest of every C source in `mexGateway.sources` when a
  build succeeds, and rebuilds only when the current digest differs -- a mechanism it adopted
  precisely because timestamps lied (a ZIP entry stores naive local time, so sources unpacked on a
  UTC+2 machine look two hours older than the binary built from them).
- On this host the digests DID differ, in both clones. The binary was genuinely stale. The warning
  was correct, and the claim that it cannot distinguish was the thing that was wrong.
- The rebuild it attempts succeeds in a cold MATLAB: 20.6 s, 257 sources, manifest updated. It
  failed during the gate runs of 2026-09-17 because the gateway was already mapped -- and those
  runs are the ones this record already documents as having had two MATLAB processes fighting over
  one DLL.

**So nothing was filed.** An issue on a dependency's tracker asserting a defect that measurement
does not support would have cost that project's owner time and this project's credibility. The
instruction was to file it; the right response to an instruction whose premise fails is to say so,
not to carry it out.

**The pattern is now four for four this month**: the accommodation models were "wrong" until
evaluated at the conditions their literature was fitted at; the `generated.tex` staleness check was
"broken" until the fault injection was verified to have injected; [H24] "gave the O-6 oracle" until
its Table 3 was read properly; and now this. **A record entry is a hypothesis until something
measures it**, and the entries most worth doubting are the ones written in a hurry at the end of a
session, which is when all four were written.

### The ninth instance, and the first where the correct value was already in the repository (2026-09-18)

O-29 recorded `Model.accel` at 88 ms per epoch (superseded) and a day at 1 Hz at 2.1 hours. **That
row had been discharged two days earlier and nobody closed it.** WP-C3 was commissioned by it, met
its claim seventyfold, and recorded the delivery in G.2 in words that were true when written:
"1.21 ms per epoch in the standard preset ... against O-29's 88 ms" -- both figures stale today.
Those two sentences sat in this file, 772 lines apart, saying different things about the same
quantity.

What that cost, in the order it happened:

- `tests/tNgfModel.m` was updated on 2026-09-15 with the true number and the true history: its
  comment records that the fast preset was 45x the cheaper until LP-C, and that the standard preset
  then cost 1.2 ms per epoch -- a figure only mildly stale now (0.798 at K=100 today).
  **The repository has therefore held the right answer, in an executed test, for three days.**
- `CLAUDE.md`, `demos/demoForceBudget.m`, `tests/tNgfApp.m` and `tests/tNgfMex.m` went on quoting
  the row. Two of them used it as a JUSTIFICATION for a design choice -- the demo's ten epochs, and
  `tNgfApp`'s use of the fast preset to avoid "18 s of gate time" that is in fact 0.17 s.
- On 2026-09-17 a user guide chapter was written fresh FROM the stale row, stating 88 ms and 2.1
  hours as the toolbox's cost, and shipped in the PDF. **A stale record does not merely fail to
  update the documents that exist; it manufactures new ones.**

And the instrument built to prevent exactly this -- `tools/benchmark.m`, written 2026-09-09 whose
own help says "a number that cannot be regenerated is a number nobody can check" -- **had never been
run, and did not work.** It errored for every value of its own argument, from the day it was
written until 2026-09-18 (see D.9).

**Then the repair nearly shipped the same defect one level down.** The first run of the fixed
instrument measured with a single warm-up call, which is not warm: it read `kernel.singleEpoch` at
2.44 ms against a settled 0.456, and the fast preset at K=1 at 31.9 ms against a settled 4.00.
Worse than the wrong digits, **an argument had already been written on top of them** -- that the two
presets differ by 4 % at K=1, so K=1 cannot be measuring physics. They differ by 4.0x. The
conclusion (batch, or use the native path) survived; the reasoning did not, and it had already been
committed to the user guide before a second run caught it. **Measuring is not one act.** A number
that has been measured once is a hypothesis about the measurement as much as about the code, and
the cheapest test of it is to measure again, differently.

**The rule this yields, which the eight earlier instances did not:** closing an obligation is part of
discharging it, and a work package that fixes a number must be the thing that deletes the old one.
Checking that a fact appears once is not enough when the second copy is an OBLIGATION to change the
first -- an open row is a standing instruction to repeat what it says. `auditPackage.py` check 13
now compares every per-epoch millisecond figure in the repository against the measured table in
`tools/benchmark_reference.csv`, so this particular fact has one authority and a machine that reads it.

### N-10a — the regression baseline 2020-01, GRACE-C and -D (issue #67, first of two PRs, 2026-10-16)

**N-10 is two PRs under #67**: **a** the frozen baseline (ORC-10: "releases change no calibrated number beyond
k x bootstrap sigma"); **b** the comparison with the hybrid ACH1B (ORC-12), which needs the ACX archives. The
classic month waits for N-08.

**Measured before scoping** (one full real day, GRACE-C 2020-01-15, from the DATAPOOL): `app.calibrate` takes
14.5 s with the `fast` preset and 86.1 s with `standard`. The fast preset's guessed exponential density gives
an along-track scale of 0.32, physically meaningless; `standard` (NRLMSISE-00) gives 0.74 and fits better (x
residual RMS 4.0e-9 against 5.1e-9). A month for both satellites in `standard` is ~1.5 h: a release-time run
on Tier A, like `tools/benchmark.m`, not a step of the 45-minute CI job.

**Scope, frozen 2026-10-16.** Tier A. Data: the DATAPOOL (`READING_DATA_POOL`; on the owner's host
`E:\DATAPOOL`), `GRACEFO\L1B\RL04\2020\2020-01-DD\`, fetched 2026-10-16 by `orbitIntegration.fetchData`
(31 archives, SHA-256 in `FETCHED.json`; the three days fetched independently on 2026-10-15 match).
- **Afterwards true:** `tools/regressionBaseline.m` (a development instrument, not a catalogued unit) calibrates
  each satellite over the month JOINTLY (N-09: one common scale, a linear bias per day, a bootstrap over days)
  against the `standard` Model along GNV1B with SCA1B attitude, and either FREEZES the result into
  `tools/baseline_reference.json` or COMPARES with it: every number may move by at most k = 3 of its sigma (the
  scale's bootstrap sigma, a bias's formal one), and the file records the toolbox version it was last verified
  at. `auditPackage.py` refuses a tree whose version differs from that record, as it refuses a stale benchmark,
  so every release must re-run the baseline. The days `run` refuses are recorded, not hidden.
- **Out of scope:** N-10b (ACH1B, ORC-12); the classic month (after N-08); running the month in CI.

- **Frozen 2026-10-16** (Tier A, R2026a Update 4, `standard`, 100 bootstrap replicates over days, k = 3;
  `tools/baseline_reference.json`):

  | | GRACE-C | GRACE-D |
  |---|---|---|
  | days used | 31 | 28 (2020-01-19..21 hold no records: safe mode) |
  | S_xx +- bootstrap sigma | 0.7620 +- 0.0084 | 0.7694 +- 0.0097 |
  | S_zz +- bootstrap sigma | 0.8248 +- 0.0019 | 0.7845 +- 0.0037 |
  | residual RMS x / y / z (m/s^2) | 4.1e-9 / 3.8e-9 / 3.4e-9 | 4.0e-9 / 1.2e-8 / 4.0e-9 |
  | run time | 2 h 9 min | 1 h 52 min |

  ~~D's cross-track residual is three times C's: its ACT1B is the GRACE-C transplant (V-5)~~ **WRONG, corrected
  2026-10-20 by N-10b:** the month's 1.2e-8 is ONE day, 2020-01-22, the first after GRACE-D's safe mode
  (19..21 Jan): alone it reads 8.5e-8 cross-track and 1.4e-8 radial, 2020-01-23 still 7e-9, and every other day
  ~4e-9, as on C. On the 26 complete days D's ACT1B cross-track rms is 4.1e-9. The transplant was a plausible
  cause written down without the per-day check that would have refuted it. Both along-track scales near 0.76 measure the `standard` Model's density
  against the instrument, not the instrument alone.
- **Found while building it.**
  1. **The first month run lost three hours**: GRACE-D's ACT1B, GNV1B and SCA1B of 2020-01-19..21 hold
     `num_records: 0`, gracefoReading refuses such a file (`gracefoReading:readL1b:noData`), and the tool wrote
     nothing until the end, so GRACE-C's finished result went with it. Now that one identifier is caught (the
     day is listed in `noDataDays`, every other reader error still stops), and each satellite is written the
     moment it is done.
  2. **Where the four hours go -- MEASURED, after first being asserted.** This record first said the cost was
     "almost all the bootstrap"; only the read (~7 min) and `calibrate` as a whole (~2 h per satellite) had been
     timed. Measured 2026-10-16 on 10 days of GRACE-C (Tier A, `standard`): read 189 s, `calibrate` 1049 s
     with `NBoot` 0 and 2049 s with 100, so **the bootstrap is 45 % of the run** (1000 s). Of the 1006 s of a
     profiled `NBoot` 0 run, **`Model.accel` takes 877 s** (1.0 ms per epoch over 864 000 epochs, in line with
     `tools/benchmark_reference.csv`), and inside it **occlusion 579 s** (`flowShadow` 489, `selfShadow` 386,
     nested) and `radiationPanelForce` 228; the CRN filter 70 s; the fit itself under 60 s. The read is
     dominated by `str2double` in the YAML parser (196 of 253 s). Two independent halves, then: the model
     evaluation (in this toolbox, the shadow lookups first) and the bootstrap (each replicate refits all days
     through `lsqvce.fit`; reusing per-day reduced normal equations would need lsqvce to hand them out or to
     accept a multiplicity per day -- its `ReducedContribution` exists internally, but `fit` exposes neither).
  3. **The comparison was fault-injected on a two-day reference** (`'Reference'` names a scratch file): an
     identical rerun is 0 sigma (bitwise reproducible), a scale moved by 4 sigma fails, other days fail.
- **Audit check 14** refuses a tree whose version is not the one the baseline was last verified at; fault-injected
  (missing file, a file verified at 0.14.0: both reported; at 0.15.0: silent).

**Delivered 2026-10-16, ☑ on Tier A**: gate **427 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 427 (no
test changed); `codeIssues` empty on the package and on the tool; `auditPackage.py` clean with check 14.

### The occlusion and radiation-force paths, made cheaper with every value bitwise unchanged (2026-10-16)

N-10a measured `Model.accel` at 877 of the 1006 s of a 10-day `calibrate`, with occlusion at 579 s. This
work package removes what can be removed WITHOUT changing a single bit of any result; what needs a
decision is listed at the end, not done.

- **Measured first** (line profile, the committed 600-epoch arc, `standard`): 1.0 ms per epoch, of which
  `geom.flowShadow` 0.49 (three shaded plates -- both inner flanks and nadir -- times 1152 quadrature
  directions, each read from the cube map at ~95 ns), `kernel.radiationPanelForce` 0.26 (577 sources x 8
  plates, each writing a strided slice of a 3xNxK array).
- **Three changes, each verified bitwise against the code on `main`:**
  1. `selfShadow`'s look-up forms the three bilinear differences once per table node instead of once per
     direction (the same subtractions of the same values; four gathers instead of five): 109 -> 88 ns per
     direction, interleaved A/B.
  2. `radiationPanelForce` accumulates one NxK array per component over all plates at once: 0.272 -> 0.195 ms
     per epoch. **A first, fully broadcast 3xNxK version was SLOWER (0.320)**: vectorising is not by itself
     faster when it multiplies the temporaries.
  3. `flowShadow` treats the eight species side by side, so its ring loop runs 72 times per block instead of
     576; that loop was the per-block cost that made a large block look profitable (block 128 -> 1024 had
     bought 21 % before it, 8 % after). The block stays 128: 8 % is not worth eight times the memory.
- **Bitwise**: `Model.accel`'s acceleration, its whole info struct and `Model.torque` on the 600-epoch arc
  equal `main`'s exactly; `flowShadow` with 'Integrals' equals it exactly; one epoch alone equals the same
  epoch in a batch. **Interleaved A/B of `Model.accel`, same host, same hour: 1.004 ms per epoch on `main`,
  0.856 and 0.841 on this branch (-15 %).**
- **The benchmark table was NOT re-committed.** Two runs on this host gave per-call costs 2-4 ms higher in
  EVERY path, including the `fast` preset, which none of the three changes touches (5.15 -> 8.97 and 7.26 ms
  at K = 1): the host was loaded. The table is check 13's authority, and a measurement of the host's load
  must not become it; it is re-measured at the next release, as every release does.
- **What remains is the look-up itself**: ~3500 directions per epoch at ~88 ns, two arctangents each. Going
  further changes numbers or principles and is the owner's call: (a) tabulate the flow-shadow ratios over
  the body-frame direction of v_rel and the speed ratio once per Model (GRACE-FO's attitude keeps v_rel in a
  narrow cone), interpolation error to be bounded; (b) evaluate `flowShadow` on a coarser time grid and
  interpolate; (c) a compiled look-up, which P11 and "MATLAB is normative" argue against.

**Delivered 2026-10-16, ☑ on Tier A**: gate **427 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 427 (no
test changed: the existing vectorisation and reference tests cover the three units); `codeIssues` empty;
`auditPackage.py` clean.

### The gas-stream shadow tabulated: geom.flowShadow Method "table", the standard preset's default (2026-10-17)

Option (a) of the occlusion work package, ordered by the owner on 2026-10-16 with "default in standard".

- **Why it can be a table.** The flow-shadow ratios depend on the geometry, the self-shadow table, the
  stream direction in the body and ONE speed ratio (every species' s is s_O sqrt(m_j/m_O)). Measured on ten
  pool days (2020-01-01/08/15/22/29, GRACE-C and -D, every 10 s, 81952 epochs): the direction stays within
  3.8 deg of +-x (C flies with -v_rel along +x, D along -x); s_O spans 7.3-10.0. GRACE-D on 2020-01-22 has
  an attitude event reaching 11 deg -- the reason the table has a guarded edge rather than a guessed one.
- **The table.** Method "lut" at the nodes of [-tan 6 deg, tan 6 deg]^2 (gnomonic v_y/|v_x|, v_z/|v_x|; 25
  nodes, 0.5 deg) x s_O in [5, 12] (15 nodes), about +x and -x; trilinear in difference form, number and
  normal clamped to [0, 1]; an epoch outside is computed by Method "lut" exactly and flagged in
  `S.tableFallback` (and so in `Model`'s info). A cone is built at first use (~4 s, 9375 nodes) and kept for
  the session, keyed by geometry, self-shadow table and quadrature.
- **Measured against Method "lut" on the 81952 real epochs:** ratios 1.2e-3 at most; drag acceleration
  6.7e-13 m/s^2 at most, rms 1.8e-13, of a 2.1e-8 drag (C; D 2.8e-13 and 9e-14); 9 of 38812 GRACE-D epochs
  fell back. Cost 2-8 us per epoch against ~450 (the per-call overhead dominates below ~100 epochs). Linear,
  not cubic: cubic halves the error but can leave [0, 1], which `dragPanelForce` refuses.
- **Found on the way:** my first prototype mirrored the table in yaw and pitch (nodes built as sg*[1; ty; tz],
  queries as v_y/|v_x|); its error was 0.1 and did NOT fall with the grid, which is what gave it away.
- **`Options.geom.flowShadow`** "exact" | "table": fast exact (no shadow at all), standard **table**, precise
  exact. `Model` refuses "table" beside `geom.selfShadow` "raycast" (the table tabulates the lut); with
  selfShadow "none" or constantCd there is no stream shadow and nothing is built.
- **The native path reads the same table**: `DynamicsSource` hands it to `mex/ngfRhs.c`, whose `ngfFlowTable`
  is the MATLAB trilinear in the same order and falls back to `ngfFlowShadow` outside it. The seven-
  configuration parity (1e-13) now covers the table through every standard-based configuration, and
  `DynamicsSourceMetamorphicFlowTable` asserts that both paths really read it (the exact shadow moves the
  result by more than 100x what the paths differ by).
- **Tests (+8):** `flowShadowTable{Contract,Reference,Precision,Vectorization,Robustness,Metamorphic,Speed}`
  and the integration test above; Options and Model refusals added to existing tests. Two of my own test
  mistakes were caught by the first run and decomposed: table-EDGE nodes come back an ulp outside the table
  (correctly computed exactly; the precision test now reads interior nodes), and a speed test drew directions
  outside the cone and timed the fallback, not the table.
- **Not related, seen on the way:** `tNgfGeometry/checkAttitudeSenseSpeed` read 17-23 against its bound 12
  on 2026-10-17 in single-class runs, on `main` exactly as on this branch, and 7.1 inside the full gate an
  hour later: the host's load, not the code (see the benchmark note of 2026-10-16). Recorded, not moved.

- **Measured after the merge, 2026-10-17** (Tier A, 10 days of GRACE-C 2020-01 from the pool, `tools/regressionBaseline`
  trial mode; `+NonGravForces` alternated between `main` before #98 (0fd75f0) and after (73c4758), one MATLAB
  process per run, order new/old/new/old then old/new):

  | `calibrate` | before #98 | after #98 |
  |---|---|---|
  | NBoot 0 | 657 s, 648 s | 385 s, 314 s |
  | NBoot 100 | 1541 s | 1204 s |

  The model half fell by 46 % (~650 -> ~350 s); the bootstrap half (~870 s) is untouched, as expected: it is
  lsqvce refitting (astrobayer/lsqvce#183). A whole run with the read (~145 s) went from 1686 to 1348 s (-20 %).
  The results move by what the table carries: S_xx 0.724226 -> 0.724230, S_zz 0.811281 -> 0.811277 (4e-6,
  0.02 of the formal sigma and 3e-4 of the bootstrap one); residual rms identical to three digits. The 2026-10-16
  split (1049 s for NBoot 0) predates #96 and was taken on a loaded host; these runs are the ones to compare.

**Delivered 2026-10-17 on Tier A**: gate **435 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 435; `codeIssues` empty; `auditPackage.py`
clean.

### The checkAlignment lags of the 0.17.0 baseline (PR #125, 2026-10-30)

The 0.17.0 baseline run raised 80 `timeSeriesPreparation:checkAlignment` warnings (ACT1B against the
CRN-filtered model inside `app.calibrate`, nothing shifted). `tools/alignmentMonth.m` (release-checks
`what: alignment`, run 41029605595, C 31 days in 2116 s, D 28 days in 1854 s) measured them per day and axis
for five variants. **Result:**

- **Not the estimator.** Removing a linear or quadratic bias per day (fitted jointly with a scale) moves no lag
  by more than 0.6 s (C) / 1.4 s (D), and the peak correlations are high: median 0.99 along-track, 0.95-0.97
  cross-track and radial. The lags are well defined; ACT1B's bias and its drift do not make them.
- **Not a time tag.** A time error is common to the axes; the measured lags are not. Along-track (drag) the
  median |lag| is 4.2 s (C) / 2.6 s (D), drifting from -12 s on 01-01 to +2..+6 s at the month's end;
  cross-track and radial they grow almost monotonically over the month, y from 3-4 s to the 60 s edge of the
  search (C from 01-26, D from 01-22), z from 4-7 s to 45 s (C) and the edge (D from 01-24).
- **In the radiation/thermal part of the model.** Against radiation + thermal alone (drag removed with the
  fitted scale), the y and z lags stay where they were (median C 17.5 / 10.5 s, D 7.6 / 17.4 s, peaks 0.95-0.97);
  the drag-only reference carries nothing on y (peak 0.06). The model's radiation + thermal signal on y and z
  is out of phase with ACT1B by seconds to more than a minute, growing over 2020-01 on both satellites.
- The fitted scales of the variants reproduce the baseline (S_xx 0.966 / 0.972, S_zz 0.820 / 0.790), and the
  warnings stay as they are: they report a real phase difference, not a defect of checkAlignment.

**Separated (PR #126, run 41042736099, C 31 days in 3888 s, D 28 in 3364 s): THE PHASE IS THERMAL INERTIA.** The
tool evaluated three models per day: the standard preset, the preset with `erp.backend = "none"` (the Sun alone)
and the preset with `thermal.model = "transient"` (Woeske's lumped walls, graceFO's heat capacities; the standard
preset's recoil is `steadyState`, without thermal inertia).

| median \|lag\| s (peak), C / D | x | y | z |
|---|---|---|---|
| standard (steadyState), bias removed | 4.2 / 2.7 (0.99) | 17.3 / 8.2 (0.97 / 0.95) | 10.5 / 18.5 (0.96 / 0.95) |
| **transient thermal, bias removed** | 7.6 / 5.0 (0.99) | **1.8 / 3.6** (0.98 / 0.96) | **2.0 / 3.4** (0.98 / 0.98) |

- With thermal inertia the y and z lags fall to 2-4 s and the correlation peaks rise; C's largest lag over
  the month drops from 60 to 10 s (y) and from 45 to 4.6 s (z). The month-long growth of the lag was the steady-
  state recoil leading the real one. A residual growth remains on D's z (0 -> 11 s over the month).
- **The calibration scale moves with it.** The per-day, per-axis scale of ACT1B on the model (one scale, linear
  bias) goes from x 0.966 / 0.972, y 0.870 / 0.903, z **0.820 / 0.790** (steadyState) to x 0.998 / 1.005,
  y 1.004 / 1.033, z **1.092 / 1.065** (transient). The S_zz ~ 0.8 of every calibration so far is, in good part,
  the steady-state thermal model, not the accelerometer.
- The along-track lag grows (-4 -> -7.6 s on C): the transient walls change x as well, and the heat capacities are
  educated guesses (graceFO's help); the x lead is the next thing to read.
- The isolated single effects were not informative for erp and thermal (peaks 0.01-0.85, lags at the edge): they
  are small and collinear with the Sun, so subtracting the others at one common scale leaves mostly the others'
  errors. srp alone carries the same lag as the whole model (it dominates y and z).
- Caveat: each day's transient model starts from the steady state at its first epoch (spin-up of the order of the
  walls' time constants is inside the day).

**The owner's rulings (2026-10-30): standard takes transient and stays integrable; the heat capacities are fitted
first.** `tools/heatCapacityFit.m` (release-checks `what: heatcapacity`, run 41140282748, C 11 days in 3923 s, D 10
days in 3425 s, every third day of 2020-01) fitted one factor in {0.25, 0.5, 1, 2, 4} on Woeske's capacity per wall
group (array, kapton, radiator) -- the lumped walls do not exchange heat, so a transient Model per factor gives every
panel's temperature and each combination is a sum of per-group recoils after the CRN filter. rms x / y / z [nm/s^2],
per-day scale + linear bias per axis, both satellites together:

| model | rms x / y / z | median scale x / y / z | median \|lag\| x / y / z s |
|---|---|---|---|
| steadyState (0.17.0) | 2.09 / 4.05 / 3.25 | 0.969 / 0.899 / 0.805 | 3.8 / 11.7 / 11.5 |
| transient, Woeske x1 | 2.21 / 3.75 / 2.28 | 0.999 / 1.021 / 1.077 | 7.0 / 2.3 / 2.9 |
| transient, best [0.5 0.25 4] | 2.16 / 3.77 / 2.20 | 0.975 / 0.961 / 0.961 | 6.7 / 2.4 / 2.1 |

- **The transient walls are the gain**: z rms falls 30 % (3.25 -> 2.28), y 7 %, x rises 6 %.
- **The capacities are NOT determined by ACT1B**: the best combination improves rms_y + rms_z by 1 % over Woeske's;
  C prefers [1 4 4], D [0.5 4 0.25] (both at grid edges); the kapton factor ranges 0.25-4 within the top five. **The
  radiator factor trades against the z scale at constant residual** (x1: S_zz 1.08, x4: 0.96), so with these data S_zz
  is known only to that band. **Woeske's values are kept** (geom.graceFO's help records why).
- **Done in PR #126 (2026-10-31):** `Options` standard sets `thermal.model = "transient"` (precise stays steadyState);
  the native path serves it through the temperature table ('ThermalOrbit'); `Model.accel` takes 'T0_K', and
  `app.forceBudget` and `app.densityInversion` carry the walls from chunk to chunk (they restarted from the steady
  state at every chunk otherwise); `NgfTestCase.standardUnscaled()` is now also steadyState (tests about something
  else); `tNgfIntegration` runs standard AS SHIPPED through 'ThermalOrbit' and compares with Model.accel along the
  whole sample orbit, and its 'Native' false comparisons take a "standardSteady" configuration. The baseline moves with
  the preset: **re-frozen 2026-10-31 (release-checks run 41156458334, 76 min, 100-replicate day bootstrap)** --

| | S_xx | S_zz | residual rms x / y / z [nm/s^2] |
|---|---|---|---|
| C, 0.17.0 (steadyState) | 0.96615 +- 0.00246 | 0.81984 +- 0.00193 | 2.13 / 3.40 / 3.46 |
| C, transient | **0.99569** +- 0.00246 | **1.09401** +- 0.00192 | 2.26 / 2.52 / 2.82 |
| D, 0.17.0 (steadyState) | 0.97193 +- 0.00267 | 0.79256 +- 0.00329 | 2.19 / 12.2 / 3.99 |
| D, transient | **1.00175** +- 0.00265 | **1.06854** +- 0.00279 | 2.27 / 12.1 / 3.09 |

  The along-track scale is now 1.00 on both satellites and the radial 1.07-1.09: 0.17.0's S_zz ~ 0.8 was the
  steady-state thermal model, and the radial scale is now bounded only by the radiator's unknown capacity (0.96-1.08 at
  constant residual, heatCapacityFit). checkAlignment still warns on every day: y and z 1-5 s, x -5 to -22 s (C) --
  the along-track lead the transient walls introduce, the next thing to read. The user guide's demo transcript (Step 4, "thermal steadyState") is stale until the
  demo is re-run.

**The x lead read at the shadow edges (PR #128, run 41187770368, standard preset as shipped = transient + conical; C 31
days in 4703 s, D 28 in 4065 s).** `tools/alignmentMonth.m` now measures the lag within +-300 s of each shadow entry and
exit (nu crossing 0.5, ~30 edges per day) and with the other thermal and eclipse model. Median signed lag in s
(positive: ACT1B lags the model), C / D, every day of the month with the same sign unless noted:

| | x | y | z |
|---|---|---|---|
| whole day (linear) | -7.6 / -5.0 | -0.7 / -3.5 | +0.8 / +3.4 |
| shadow ENTRY, conical | **-3.9 / -4.4** (31/31, 28/28 negative) | -4.7 / -6.7 | -3.9 / -3.7 |
| shadow EXIT, conical | **+1.0 / +1.2** (all positive) | +0.8 / +0.6 | +3.2 / +3.7 |
| entry, solaarsCF | +1.1 / +0.4 (mixed sign) | +0.9 / -1.9 | +0.8 / +1.0 |
| exit, solaarsCF | **-8.7 / -8.7** (all negative) | -9.4 / -10.3 | -7.0 / -6.8 |
| entry / exit, steadyState walls | -3.7 / +1.1 (C) | -3.9 / +0.6 | -4.0 / +4.9 |

- **It is the shadow, and it is not a time tag.** ACT1B sees the Sun go ~4 s EARLIER at entry and come back ~1-3 s LATER
  at exit than the conical model says -- on all three axes and on every day. A time-tag error would shift both edges the
  same way. The real shadow is WIDER than the geometric cone (~4 s x 7.6 km/s ~ 30 km at entry), which is what the
  atmosphere's extinction and refraction do; the thermal model does not move the edges (steadyState gives the same).
- **solaarsCF fixes the entry and breaks the exit**: entry ~0 s, but exit -8.7 s on every day of both satellites, i.e.
  solaarsCF holds the shadow ~9 s too LONG after exit. An asymmetry of that size between entry and exit, uniform over a
  month, is more likely an implementation question than physics (the attenuation curve is a function of the grazing
  height, which entry and exit traverse in opposite directions): next is a symmetry check of `env.eclipseFactor`
  solaarsCF on mirrored geometry against `tools/pyxcheck`'s coefficients.
- The whole-day x lag (-7.6 s) is not the edges alone: with steadyState walls it is -4.2 s; the transient walls add
  ~3 s of their own on x (kapton front/rear), which the heat-capacity fit already showed.
- **precise (run 41187772070) stopped after three days of GRACE-C**: the committed CERES excerpt ends at 2020-01-16 and
  `earthRadiance` refused the next day (`release-checks` now runs precise on 1:3:15, and the tool refuses later days up
  front). The three days (01, 06, 11; precise = steadyState + solaarsCF + unscaled msis21) say the same as standard's
  edges: entry -10 / -4 / -1 s, exit -12 / -8 / -9 s on y/z -- solaarsCF's exit again ~9 s late; whole-day y/z
  -0.2/+0.6, +5.7/+5.0, +10.6/+7.9 s steady against -8.9/-9.0, -5.4/-5.6, -2.6/-3.3 s transient. With solaarsCF's exit
  error in the same signal the thermal choice for precise cannot be read cleanly: **solaarsCF first, then precise.**

**solaarsCF's oblateness mapping was INVERTED (PR #129, 2026-10-31).** `env.eclipseFactor`, `mex/ngfRhs.c` and the
oracle's docstring stretched the IN-ECLIPTIC component of rE by sO = a/b, i.e. modelled an Earth flattened at the
equator. The silhouette's short semi-axis lies along the ecliptic normal (at the January solstice exactly the projected
spin axis), so the component ALONG the normal is the one to stretch onto the fit's sphere of radius a. A stdlib
simulation (circular 490 km orbit, i = 89 deg, January pole) gives solaarsCF minus conical at the 50 % crossing: inverted
mapping entry -8.0 / exit +8.0 s at beta 0, -8.3 / +8.6 s at beta 20, -10.9 / +5.4 s at beta 40 -- the size and the
asymmetry ACT1B measured (+5.0 / -9.7 s against conical); corrected mapping -0.4 / +0.4 s, -0.8 / +0.4 s, -0.7 / +6.2 s.
The physical anchor, now a test (`tNgfEnv/eclipseFactorMetamorphic`): the half-shadow limb sits the same 22.5 km above
the equatorial radius in the ecliptic plane and above the polar radius along the normal (the inverted mapping put the
polar limb 42.5 km above b). The thesis itself (VTechWorks) was not reachable from the session; the decision rests on
the geometry and on the measured asymmetry. Affects only srp.eclipse "solaarsCF" (precise; standard is conical).
Predicted gate unchanged at 511 (a new bound inside an existing method). Next: re-run alignment for standard
(altEclipse = solaarsCF) and precise (1:3:15) to read the corrected edges, then decide precise's thermal model.

**The corrected solaarsCF against ACT1B (run 41213058839, standard + altEclipse solaarsCF with the PR #129 mapping;
C 31 days, D 28).** Median signed lag at the shadow edges, s (C / D):

| | x entry | x exit | y entry | y exit | z entry | z exit |
|---|---|---|---|---|---|---|
| conical (the preset) | -3.85 / -4.38 | +1.00 / +1.24 | -4.70 / -6.67 | +0.81 / +0.60 | -3.93 / -3.68 | +3.23 / +3.72 |
| solaarsCF, inverted (run 41187770368) | +1.13 / +0.39 | -8.67 / -8.68 | +0.86 / -1.93 | -9.44 / -10.30 | +0.80 / +1.04 | -7.02 / -6.75 |
| **solaarsCF, corrected** | **+1.09 / +1.31** | **+0.49 / +0.64** | **+0.37 / -0.26** | **+0.83 / +0.67** | **+0.96 / +1.32** | +2.98 / +3.63 |

- The corrected SOLAARS-CF puts both edges within ~1 s of ACT1B on x and y (largest |lag| at entry 1.9 / 2.0 s on x
  against conical's 11.8 / 11.6; days aligned within 1 s 13 and 22 of 31 on x for C against 3 and 16). z keeps +3 s at
  EXIT under every eclipse model and both thermal models' windows -- that is the radial thermal response after sunrise,
  not the shadow. **The atmosphere's widening of the shadow is real and SOLAARS-CF, correctly mapped, carries it.**
- Whole-day lags move little and not uniformly (C y -0.7 -> +2.2, z +0.8 -> +2.8; D y -3.5 -> -0.1, z +3.4 -> +6.6; x
  -7.6 -> -7.0 / -5.0 -> -4.2): over a whole day the thermal response, not the edge timing, sets the phase.
- **Owner's call:** whether `standard` takes `srp.eclipse = "solaarsCF"` (a preset change: native path supports it,
  baseline re-freeze) -- the edges argue for it, the whole-day lags neither for nor against.

**precise with the corrected solaarsCF (run 41213060789, every third day 2020-01-01..13, C and D; precise = steadyState
walls + solaarsCF + unscaled NRLMSIS 2.1).** Median signed lag, s (C / D):

| | y whole day | z whole day | x entry | x exit | z entry | z exit |
|---|---|---|---|---|---|---|
| precise as shipped (steadyState) | +13.7 / +11.6 | +12.0 / +19.5 | +0.4 / +0.6 | +1.1 / +1.2 | -1.3 / -0.2 | +5.9 / +7.1 |
| precise with transient walls | **+3.4 / +1.5** | **+2.6 / +6.7** | +0.3 / +0.5 | +0.9 / +1.1 | -0.8 / +0.1 | **+3.1 / +3.7** |
| precise with conical | +8.2 / +5.7 | +7.1 / +13.1 | -7.3 / -7.1 | +1.7 / +1.7 | -8.9 / -8.3 | +5.2 / +6.0 |

- **The steady state leads ACT1B by 12-20 s on y/z in precise too; transient walls take it to 2-7 s**, as in standard.
  The x whole-day lag sits at the -60 s edge: precise's density is unscaled (the factor of 2 in solar minimum), so the
  drag phase is not what is measured there. With solaarsCF the edges are within ~1 s; conical's entry is 7-9 s late.
- **Recommendation for the owner: precise takes thermal.model "transient"** (RULED 2026-10-31: yes; standard keeps conical) (no baseline is frozen for precise; the native
  path serves it through 'ThermalOrbit').
- CI on a086fa7 was red BEFORE the gate: Code Analyzer ISCL in the new `igrfCoefficients` refusal (`numel(...) == 1`
  instead of `isscalar`) -- fixed in the next commit; `codeIssues` cannot run here, so a new line only shows on CI.
  CI on fa9707a: 510/511, `tNgfAttitude/attitudeProviderItsgContract` errored -- `[tc.f.itsg tc.f.itsg]` joined two
  char vectors into ONE nonexistent path (mustBeFile); now a two-element string array. ☑ CI 41245260551 on fe5c8b3:
  511/511, predicted 511.

**Next, the owner's call:** whether `standard` should take `thermal.model = "transient"` (it changes the baseline's
S_zz from ~0.8 to ~1.08 and needs a re-freeze), and whether the heat capacities are to be fitted (the lag per
panel group is the observable).

### GRACE-D 2020-01-22, the day that is D's cross-track residual (PR #131, 2026-11-01)

The frozen baseline's cross-track residual of D over 2020-01 is 12.08 nm/s^2 against C's 2.51. The per-day bias sigmas of
`tools/baseline_reference.json` (the fit estimates the noise per day) put it on one day: 2020-01-22, the first day after
the safe mode (19-21 have no ACT1B), ~85 nm/s^2 over about half a day; without it the month is ~3.2 nm/s^2. Every other
day of D is 2.6-3.4 nm/s^2 in `dragScaleMonth` (run 40849007121, firings cut), which skipped 01-22 for its non-GPS THR1B.

`tools/dayDiagnosis.m` (run 41325337776; standard preset unscaled, one bias per day, 10 s), per GPS hour:
- **D 01-22:** ACT1B 43015 records from ~12:03. Hour 12: rms y 293.8 nm/s^2, mean y -21.6, ACT1B 1 Hz noise y 14.0,
  112 firings / 90.6 s on-time, SCA1B rate scatter 192.6 urad/s. Hours 13-23: rms y 3.9-12.9, noise 0.35-0.97, 5-25
  firings / 0.5-2.5 s, rate 7.4-9.1. Day: rms y 83.18; **with every firing cut 0.75** (26 % of epochs cut). THR1B time_ref
  G:283 R:27; the 27 R records are in hours 6-11, before ACT1B resumes.
- **D 01-23:** rms y 9.42 (hours 0-13: 7.0-21.3, noise 0.6-1.6, 12-44 firings; hours 14-23: 3.0-5.4, noise 0.3-0.4);
  firings cut 0.79 (38.5 % cut). 141 flagged ACT1B records in hour 17. Hour 10: z 15.05 with rate 24.0 urad/s.
- **C 01-22 (reference):** rms y 2.46, noise 0.20, 6-23 firings / 0.3-1.2 s per hour, rate 2.9-5.3; firings cut 0.70.

Reading: the excess is attitude re-acquisition on thrusters after the safe mode, concentrated in the first hour of data
and decaying over ~26 h; outside the firings ACT1B agrees with the model as well as C does. The THR1B time_ref is not the
problem (the R records precede the data). **`tools/regressionBaseline.m` calls `app.calibrate` without 'Thrusters'**, so
the frozen baseline fits through every firing on both satellites; D simply fires more. Next (owner's call): pass THR1B to
the baseline's calibration (readTHR1B must then accept a file whose non-GPS records lie outside the data, or drop them),
re-freeze; or exclude the first hours after an ACT1B gap.

### Release 0.18.0 (2026-10-31, PR #130)

What it carries: the transient thermal walls in `standard` (PR #126) and `precise` (PR #129, owner), served on the native
path by a temperature table ('ThermalOrbit'); SOLAARS-CF's oblateness mapping corrected (PR #129); transient walls carried
across `Model.series` / `forceBudget` / `densityInversion` chunks (PR #128); the eight defects of the help audit (PR #129,
among them the exact FNV-1a in `provenance`, so provenance hashes differ from 0.17.0's); the generated test catalogue and
the documentation pass (PR #128). `standard` keeps `srp.eclipse = "conical"` (owner, 2026-10-31).

Release measurements on the CI runner (`release-checks.yml`, head 63e6a63):
- **Baseline PASSED** (run 41248944129, 76 min): every compared number of 2020-01 C and D (usedDays, S_xx, S_zz, the
  worst of 186 / 168 biases) is 0 sigma from the freeze of 2026-10-31 01:42 -- the chunk-overlap fix and the other
  changes do not move the calibration at the printed digits. `meta.toolbox` 0.18.0, `meta.verified` 2026-11-01 00:11 UTC.
  checkAlignment warned 60 times (31 C, 29 D; C lags up to ~10 s, D x -2..-17.5 s, one D y at the 60 s edge).
- **Benchmark** (run 41248946430): `Model.accel` standard 26.5 / 3.99 / 0.762 ms per epoch at K = 1 / 10 / 100 (0.17.0
  on the same runner: 20.0 / 3.17 / 0.650 -- the transient walls' sub-steps, +17 % at K = 100), fast 14.7 / 1.51 / 0.247,
  `kernel.singleEpoch` 1.30 ms (0.895), native `DynamicsSource.accel` 1.31 ms (0.996). The compiled paths did not change
  in code for the standard configuration measured; their +40 % is the runner's own spread (cf. fast 0.119 -> 0.241 at
  0.17.0 with no code change) and is recorded, not explained.
- Gate 511/511 (CI 41247178221 on the merged #129 head).

### Release 0.17.0 (2026-10-29)

MINOR (added functionality, VERSIONING.md) with one deliberate preset change, by owner ruling: the `standard` preset
takes `drag.densityScale = "combined"` (PR #122), strictly -- a standard Model evaluates only with 'DensityScale'. Also
released: DTM2020 (opt-in, CNES licence), `io.readSwarmDensity` (Swarm DNSxPOD and ESA TOLEOS), `env.densityScale`,
`env.densityDriverScale`, the scale on the native integrator path, Dst/F30 reading, the pool fetches (#107-#122).
Every version copy bumped. **The two release measurements ran on the CI runner for the first time**
(`.github/workflows/release-checks.yml`, dispatched on the release branch), not on Tier A:

- **Regression baseline RE-FROZEN** (run 40911119628, Mode "freeze", 75 min on the runner against 3.1 h on Tier A
  at 0.16.0), because the preset change moves the calibration on purpose -- a "compare" against 0.16.0's numbers
  would fail by design. Joint over the month, 100-replicate day bootstrap, standard preset with the combined
  scale from the pool's Swarm month and CLS F30:
  C (31 days) S_xx **0.96615** +- 0.00246, S_zz 0.81984 +- 0.00193 (0.16.0: 0.762006, 0.824845);
  D (28 days, 19-21 without records) S_xx **0.97193** +- 0.00267, S_zz 0.79256 +- 0.00329 (0.16.0: 0.769365, 0.784481).
  **The along-track scale moves from 0.76-0.77 to 0.97 on both satellites**: with the model drag halved by the
  observed density, ACT1B no longer has to absorb a factor of 2 in its scale -- the 0.76 of 0.16.0 was the MSIS
  density error, not the accelerometer. Radial S_zz moves by under 1 % (C) and +1 % (D), the cross-track is held.
  Residual RMS x / y / z: C 2.13 / 3.40 / 3.46 nm/s^2, D 2.19 / 12.2 / 3.99 (D's y carries 2020-01-22, N-10a).
- **Benchmark** (run 40914918235; the first run 40911108153 failed: the workflow did not build the kernel MEX it
  times, fixed in b504ae2): standard `Model.accel` 20.0 ms at K = 1 and **0.650 ms at K = 100**; fast 13.3 and 0.241;
  compiled kernel 0.895 ms; native `DynamicsSource.accel` 0.996 ms. **Host change, not only code change**: 0.16.0
  measured 11.9 / 0.440 on the owner's host, and the fast preset -- untouched by this release -- moved 0.119 -> 0.241 at
  K = 100 (x2.0) while standard moved x1.5, so the runner is the whole difference and the density scale costs nothing
  visible. The user guide and CLAUDE.md quote the new table.

### Release 0.16.0 (2026-10-18)

MINOR (added functionality, VERSIONING.md): N-05 figures, N-09 multi-arc calibration, N-11a regressors, the N-10a
baseline instrument and audit check 14, the bitwise-neutral speed-up of #96 and the flow-shadow table of #98.
Every version copy bumped (toolboxVersion, the seven Contents.m, README, CLAUDE.md, HANDOVER's title, the
generated doc pages). Issues #59, #60, #62 and #66 closed with delivery notes the same day (their PRs had
carried no closing keyword).

- **Gate** (Tier A, fresh dependency clones, R2026a Update 4): **435 passed, 0 failed, 0 incomplete, GATE TRUE**,
  predicted 435.
- **Benchmark** (`tools/benchmark.m`, quiet host, committed `tools/benchmark_reference.csv`): standard
  `Model.accel` 11.9 ms at K = 1 and **0.440 ms at K = 100** (0.15.0: 17.0 and 0.836); fast 5.64 and 0.119;
  compiled kernel 0.487 ms; native `DynamicsSource.accel` 0.596 ms (0.505 at 0.15.0 -- one run, not decomposed;
  the parity tests bound its values, not its time). The user guide's cost section and CLAUDE.md quote these.
- **Regression baseline, the first release it gates** (`tools/regressionBaseline`, Mode "compare", the whole
  month from the pool): **pass**. Every quantity within 3 sigma, the worst 0.013 sigma (a GRACE-D bias); scales
  moved by at most 0.0012 sigma (C S_xx 0.762006 -> 0.762008, S_zz 0.824845 -> 0.824843; D S_xx 0.769363 ->
  0.769365, S_zz 0.784481 -> 0.784480) -- the flow-shadow table's footprint, as its measurement predicted.
  The same days were used and refused. `tools/baseline_reference.json` now carries `toolbox 0.16.0` and the
  verification time; its frozen numbers are untouched. The month took 11193 s (3.1 h) against 14454 s at
  0.15.0: C calibrate 7278 -> 5573 s, D 6318 -> 4733 s; the bootstrap is now most of it (lsqvce#183).
- `auditPackage.py` clean after the gate, the benchmark and the baseline; `codeIssues` empty. O-35 closed.

### The day bootstrap on lsqvce's per-day contributions (lsqvce#183, 2026-10-19)

lsqvce 6.35.0 (its WP-158, PR astrobayer/lsqvce#210, closing #183, requested by this toolbox on 2026-10-17) hands
out every day's reduced contribution of a fit (`R.contributions()`) and solves any weighted combination of days
(`lsqvce.solveContributions(C, c)`), one replicate in O(D p^2). `app.fitCalibration`'s common mode now fits once,
takes the contributions, and makes each replicate the drawn days' multiplicities on them (`localReplicate`); the
refit per replicate is gone. Minimum lsqvce raised from 6.22 to 6.35.

- **Under the defaults it is the same estimator.** White noise and no robust weights: the replicate IS the joint
  refit of the drawn days. Measured on 10 real days of GRACE-C (2026-10-19, `calibrate` with 100 replicates,
  alternating with `main` new/old/new): S_xx +- 0.01386824656 and S_zz +- 0.002115559478 on both, identical in
  all ten printed digits; `calibrate` **1410 s -> 329 s and 356 s**, i.e. the bootstrap went from ~1060 s to
  noise beside the ~350 s fit. The test `fitCalibrationMultiArcDayBootstrapPrecision` checks it against an
  independent oracle: for the very draws `timeSeriesPreparation.bootstrap` makes, the public joint fit over the
  EXPANDED day list (a day drawn twice generated twice at another day's time; `syntheticPair` is shift-invariant)
  -- 2.2e-12 relative in sigma. `fitCalibrationDayBootstrapSpeed`: 50 replicates cost 1.03x the fit without them.
- **Under 'Noise' "estimate" or 'Robust' it is not.** The variance components and the robust weights are now held
  at the full fit's values in every replicate (lsqvce: fixed weights); before, every replicate re-estimated them.
  A bootstrap under fixed weights is a standard, defensible estimator, but it is a change -- Part I, open to the
  owner. The single-arc block bootstrap still refits, since it resamples samples, not days.

**Delivered 2026-10-19 on Tier A**: gate **437 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 437 (+2);
`codeIssues` empty; `auditPackage.py` clean.

### N-11b — the temperature-dependent scale (issue #68, second of two PRs, 2026-10-20)

CONCEPT 6.2: "an optional temperature-dependent, fully-populated scale matrix". It builds on N-11a's regressors.

- **Afterwards true.**
  - `app.fitCalibration` takes `'ScaleTemperature'`, the name of one `'Regressors'` column T. Every ESTIMATED
    scale entry (the 'Scale' mask: diagonal, or fully populated with "full") gains dS_ij/dT, from the column
    H((T - T0) a_j): the series linearly onto the raw grid, times the model axis, then the filter.
  - S(T) = S + S1 (T - T0), where `scale` is S at T0. T0 is the mean over the raw grid of the arc, and over
    ALL arcs in common mode, where S1 is common like S.
  - The new columns enter the observability gate. A VIF above 'VifMax' gives "warn", not a refusal.
  - It is refused beside a time-varying scale (`unsupportedScaleTemperature`), and a name that is not a
    regressor column is refused too (`inputScaleTemperature`).
  - `app.Calibration` carries `scaleTemperature` (.name .units .reference .coef 3x3 .cov 9x9), schema "4"
    (schemas "1" to "3" load).
    - `apply` evaluates S(T(t)) per epoch and refuses without the series.
    - `invert` refuses (S(T)^-1 is not linear in T); so does a calibration combining it with a scale curve.
    - New methods: `isTemperatureDependent`, `scaleAtTemperature`, `temperatureAt`; `summaryTable` gains
      dS_ij_dT rows.
  - `app.syntheticPair` plants 'ScaleTemperature' on the first regressor column, optionally with a fixed
    'ScaleTemperatureReference'.
- **Why one reference over all days.** GRACE-FO's accelerometer is thermally controlled. TSU_Y+ moved 0.065 C
  over ten days and its daily means 0.03 C (measured below). A per-day reference would have removed exactly
  the day-to-day differences that are S1's only leverage across days. The synthetic oracle shows the common
  reference working on daily offsets alone: three days whose in-day swing is 5 mK return S1 to 2.4e-11
  relative.
- **Found while building it.**
  1. **N-11a's interpolation bridged gaps.** `Calibration.interpolateSeries` interpolated linearly across
     ANY gap in a measured series. A 20-minute housekeeping gap moved a planted S1 by 3.7 %, and the rows
     inside the gap were still fitted, on invented temperatures. It now returns NaN where the samples around
     an epoch lie more than 10x the 90th percentile of the series' sample spacing apart. This is one rule for
     the estimator and for `apply`, and the regressor terms of N-11a get the same fix.
  2. **The first version of that rule used the MEDIAN spacing, and the real-data test refused it.** Real
     series come in bursts, so their median spacing is the burst's:
     - MAG1B: pairs 0.5 ms apart, once a second (median 0.5 ms, 90th percentile 1.0 s);
     - AHK1B's TSU: 0.1 s bursts every 4.7 s.
     Every epoch between bursts became a "gap", and TSU_Y was refused as having no finite value.
  3. **`apply` with S(T(t)) does not commute with the filter** the fit put around S a. The difference is
     about |S1| times the temperature change within one 140.7 s window: 2.2e-4 relative for a deliberately
     exaggerated |S1| = 2e-3 per C on a 0.5 C orbital swing. On GRACE-FO (0.02 C a day) it is negligible.
     The help names and quantifies it, and the test bounds it by that formula.
- **Measured on real data** (Tier A, 2026-10-20). Ten days of GRACE-C, 2020-01-01..10, from the pool;
  `calibrate` with the standard preset, 100 replicates, TSU_Y+ from AHK1B as regressor and as the scale's
  temperature.
  - dS_xx/dT = 0.26 +- 0.83 per C (0.3 sigma); dS_zz/dT = -0.13 +- 0.13 per C (1.0 sigma).
  - Gates pass (VIF 1.2); residual rms unchanged to four digits.
  - S_xx 0.7248 +- 0.0138 without, 0.7246 +- 0.0146 with.
  - **On thermally controlled GRACE-FO a temperature-dependent scale is not observable in ten days**, and the
    toolbox now says so with a number rather than an assumption. Its use is a mission or epoch with real
    thermal excursions (classic GRACE after its thermal control changed, V-10 / N-08).
- **Tests (+6):** `fitCalibrationScaleTemperature{Reference,Contract,Robustness,Metamorphic}` and
  `calibrationScaleTemperature{Contract,Precision}`. `syntheticPairContract` and `calibrationContract` were
  updated for the new truth fields and schema "4".
  - noise-free S1: 2.8e-13 (one arc), 2.4e-11 (three days, offsets only);
  - temperature in another unit: S1/10, scale unchanged;
  - across a 20-minute temperature gap: rows dropped, S1 returned;
  - `apply` against the per-epoch observation equation;
  - `apply` of a fitted calibration within 2 |S1| dT_window.

**Delivered 2026-10-20, ☑ on Tier A**: gate **443 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 443
(437 + 6); `codeIssues` empty; `auditPackage.py` clean. **N-11 (#68) is complete.**

### N-10b — the hybrid transplant ACH1B as an oracle for GRACE-D (issue #67, second of two PRs, ORC-12, 2026-10-20)

ORC-12: "calibrated ACT1B compared with the hybrid ACH1B (ACX/ACX2) for the same days". ACH1B (from the `.ACX`
daily archives) exists for GRACE-D only. It is ACU + ACS + ACM + AC0 (CONCEPT V-16): the updated transplant from
GRACE-C plus the SDS's own corrections for D's radiation pressure (ACS), mass transfer (ACM) and thrusts (AC0).
It is linear accelerations only.

- **Data.** Approved by the owner: the 31 `.ACX` archives of 2020-01 from GFZ ISDC, 460 MB, into the DATAPOOL.
  Each archive's SHA-256 is in `2020/FETCHED.json`, and every extracted file was verified against the archive's
  own MD5 `checksum`. 19..21 Jan are empty archives (the safe mode), and 18 and 22 Jan are partial.
- **Afterwards true.**
  - `tools/regressionBaseline` takes 'Product' "ACH1B" (GRACE-D, trial mode only; only ACT1B is frozen and
    compared). Its arcs are LABELLED "ACT1B": timeSeriesPreparation's `Config` has no ACH1B preset, and the
    hybrid product is CRN-filtered and flagged like ACT1B, so the label prepares both alike.
  - `io.readACC1B`'s help no longer calls ACH1B "high rate".
  - 120 real records of each of the five `.ACX` products are fixtures in `data/samples/acx`, and
    `readACC1BHybridReference` asserts ACH1B = ACU + ACS + ACM + AC0 through gracefoReading's independent generic
    reader: 7.5e-16 against a computed 2.9e-15, in a window where the thruster model is non-zero.
- **Measured** (Tier A, standard preset, 100 day-bootstrap replicates; the 26 days on which both products are
  complete; the same code for both):

  | GRACE-D | ACT1B | ACH1B | |
  |---|---|---|---|
  | S_xx | 0.7688 +- 0.0096 | 0.7721 +- 0.0098 | 0.24 sigma apart |
  | S_zz | 0.7841 +- 0.0037 | **0.8147 +- 0.0026** | **6.8 sigma apart** |
  | residual rms x / y / z | 3.96 / 4.12 / 3.61 e-9 | 4.16 / 4.12 / 3.07 e-9 | z -15 %, x +5 % |
  | bias, median day | | | z -1.5e-8, x -1.6e-9 m/s^2 |

  - **Radial.** The hybrid product moves D's radial scale from 0.784 to 0.815, towards C's 0.825 (the frozen
    baseline), and lowers the radial residual by 15 %. Against this toolbox's independent force model, the
    SDS's radiation-pressure correction (ACS) is doing what it is meant to do.
  - **Along-track and cross-track** agree between the products within their sigma.
  - The along-track residual is 5 % higher with ACH1B, which is not explained yet.
- **Found while building it.**
  1. **The N-10a attribution of D's cross-track residual to the transplant was wrong** (corrected in place
     above): it was 2020-01-22, the day after the safe mode. A monthly rms can be one day; per-day numbers
     first, then a cause.
  2. **`regressionBaseline` failed on a single day**: `app.calibrate` returns its single-arc result, without
     `used`, for one arc. It is fixed in the tool.
  3. **timeSeriesPreparation has no ACH1B preset**: `Config` maps ACC1B/ACT1B to the ACC CRN preset and mask
     bits [5 6 7]; ACH1B gets neither. It is worked around here by the label, and belongs upstream.
     **Closed upstream the same day**: timeSeriesPreparation 1.5.0 (astrobayer/timeSeriesPreparation#74, its
     X-19) prepares ACH1B exactly as ACT1B. The label is gone from `regressionBaseline`, and the minimum is
     raised from 1.4 to 1.5. The label was equivalent: 2020-01-05/06 of GRACE-D on ACH1B gave scale, bias,
     bias sigma and residual rms **bitwise** identical with the label and with the product's own name, so
     the measurement above stands. Gate 444/444, predicted 444.

**Delivered 2026-10-20, ☑ on Tier A**: gate **444 passed, 0 failed, 0 incomplete, GATE TRUE**, predicted 444
(443 + 1); `codeIssues` empty on the package and on the tool; `auditPackage.py` clean. **N-10 (#67) is complete.**

### N-11a — measured regressors in the calibration: temperatures and magnetorquer currents (issue #68, first of two PRs, 2026-10-15)

**N-11 is two PRs under #68**: **a** a general 'Regressors' term (CONCEPT 6.2: "a temperature regressor in the
bias basis (housekeeping product per V-13) ... and magnetic-torquer current as a regressor"); **b** the optional
temperature-dependent, fully-populated scale matrix, which builds on a.

**Measured before scoping (real AHK1B, 2020-01-01 GRACE-C, gracefoReading's generic reader).** The GRACE-FO
accelerometer is thermally controlled: over the whole day the sensor-unit temperatures TSU Y+/Y-/Z+/Z- span
0.022 C (sigma 0.006), the front-end electronics TFEEU 0.04 C, the control unit TICUN 0.85 C; TICUR reads
-273.20 on every record (a dead channel). Each temperature arrives about every 2.4 s, with `prod_flag` choosing the
columns per record. So within ONE day a temperature regressor has almost no leverage against the bias; it is a
campaign quantity (N-09), and the help says so. The magnetorquer currents of MAG1B, by contrast, swing by
amperes every orbit.

**Scope of N-11a, frozen 2026-10-15.** Tier A.
- **Afterwards true:** `app.fitCalibration` takes 'Regressors': a timeSeriesPreparation.Arc of K measured series
  on any grid (1xN for several arcs, or one covering all). Each is interpolated LINEARLY onto the kernel's raw
  grid (no overshoot between housekeeping samples), demeaned over the arc's raw span, and filtered by H like the
  bias basis; it enters every axis. Single arc: its columns join the bias block. Several arcs in common mode: its
  coefficients are COMMON (a physical coupling does not change from day to day). A series with no finite value,
  or constant over the span, is refused (`inputRegressors`); TICUR is NOT caught by that (it wobbles 3 mK about
  absolute zero), so choosing live channels is the caller's, and the help says so. `app.Calibration` carries `regressors`
  (.names .units .mean .coef 3xK .cov) and schema "3"; `apply(arc, 'Regressors', R)` subtracts sum_k c_ik
  (z_k(t) - mean_k) and REFUSES a calibration with regressor terms when R is not given; `invert` maps the
  coefficients as it maps the bias. `app.syntheticPair` plants a coupling ('RegressorFcn', 'RegressorCoupling')
  and returns the series as `truth.regressors`.
- **Out of scope:** N-11b; choosing a channel for the user; a nonlinear (e.g. quadratic) temperature term.
- **Oracles.** A noise-free pair with a planted coupling returns it to the computed rounding bound; a regressor
  scaled by k has its coefficient divided by k; apply with the regressor series undoes the planted coupling;
  real data: 30 min of 2020-01-01 GRACE-C with TSU_Y+ and the MAG1B torquer currents calibrates end to end.

- **Found while building it.**
  1. **TICUR is not caught as dead by a constancy check**: on 2020-01-01 it reads -273.205 to -273.202, absolute
     zero plus a 3 mK quantisation wobble. The scope claimed it would be refused; it is not, the help says
     choosing live channels is the caller's, and the real-data test records TICUR's distance from absolute zero.
  2. **Provenance JSON-encodes the options**, so a function handle (syntheticPair's 'RegressorFcn') and an Arc
     (fitCalibration's 'Regressors') each need a hashable stand-in: the function's text, the regressor names.
  3. **The first gate was 426/427**: `syntheticPairContract` pins `truth`'s field list, and the two new fields were
     not added to it. Fixed in the test (with an assertion that nothing is planted by default); re-run 427/427.
- **Measured on real data** (30 min, 2020-01-01 GRACE-C): the fit with TSU_Y+ and the three net magnetorquer
  currents completes with finite coefficients; TSU_Y+ moved < 0.05 C in those 30 min.

**Delivered 2026-10-15, ☑ on Tier A** (`ciGate` over fresh clones, R2026a Update 4): **427 passed, 0 failed,
0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly: 421 + 6 (`tNgfCalibration` regressor contract and
precision; `tNgfFitCalibration` regressor reference, robustness, metamorphic; `tNgfCalibrate` real day).
`codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when re-run after the gate. 59 units.

### N-09 — multi-arc calibration: one scale over many days (issue #66, 2026-10-15)

**Scope, frozen 2026-10-15.** Tier A. CONCEPT 6.6 ("`app.calibrate` accepts an `Arc` array: a multi-arc
calibration. Per-arc biases are lsqvce local parameters, pre-eliminated; the scale is common or per arc as
configured ... Multi-arc fits use lsqvce's streaming accumulation (declaration kind A, one system per arc) ...
and keep the three axes as three separate problems"); the v3 plan's WP-19 part.
- **Afterwards true:** `app.fitCalibration(obs, ref)` takes 1xN Arc arrays (one per day; `ref` 1xN or one arc
  covering all). 'ScaleMode' "common" (default for N > 1): per axis ONE lsqvce kind-A system whose days are the
  arcs and whose arcs are each day's gap-free runs; the scale (and K2) columns are COMMON (A), each day's bias
  basis is a 'day'-scope LOCAL block (B), pre-eliminated by lsqvce and recovered by back-substitution
  (`Results.residuals().sys.y`, `.Cyy`). "perArc": N independent single-arc fits. It returns a 1xN
  `app.Calibration` array, each with the common scale and its own day's bias. Common mode refuses what it would
  only half-support: a knotted bias basis (a day's local block must have one size), 'Jumps', and a
  time-varying common scale (`unsupportedMultiArc`). Uncertainty: the common scale's covariance is a bootstrap
  over DAYS (`timeSeriesPreparation.bootstrap` on the day index, block length 1, each replicate refitted
  jointly with its days relabelled): days are independent where the samples inside a day are not. It warns
  below 10 days. The per-day biases carry lsqvce's formal local covariance, and the help says so.
  `app.calibrate(acc, ...)` with an Arc array runs `timeSeriesPreparation.run` in campaign mode, builds each
  day's reference (Model along each day's orbit, or a reference Arc per day), fits jointly, and leaves out a
  day run refused, recording why.
- **Out of scope:** N-10 (the month baseline); a common time-varying scale; per-day knotted biases and jumps in
  the joint fit; a weighted bootstrap.
- **Oracles.** "perArc" equals N single-arc fits, bitwise; a noise-free multi-day pair with ONE planted scale
  and a different bias per day returns them to the computed rounding bound; the joint fit with the bias
  eliminated equals the dense joint least squares (a block-diagonal bias design next to the common scale
  columns) solved by QR here in the test (the algebra of pre-elimination, at rounding); real data: three days of
  GRACE-C 2020-01-01..03 (30 min excerpts of ACT1B, GNV1B, SCA1B, fetched from GFZ ISDC on 2026-10-15; archives,
  sizes and SHA-256 in `data/samples/README.md`) calibrate end to end with one common scale.
- **Found while building it.**
  1. lsqvce wants each arc's `B` to span EVERY local column of the system (m_a x nl over all days), the arc's
     own day block in place; the docs say so and a first attempt that passed only the day's block was refused.
  2. An axis whose scale is held (S_yy by default) has NO common column, and lsqvce takes no system without
     one: those days decouple into independent bias fits, which go through the single-arc solver.
  3. Measured on real data (30 min of each of three days, the fast Model): every day passed `run`, the common
     scale is one matrix in all three Calibrations, and the x biases of consecutive days agree to 8 %
     (-3.38, -3.32, -3.12e-7 m/s^2).
- `fitCalibration`'s vectorization exemption ("a campaign of arcs is a different estimator (N-09)") is DELETED:
  "perArc" mode equals N single-arc fits, bitwise, and that is now its vectorization test.

**Delivered 2026-10-15, ☑ on Tier A** (`ciGate` over fresh clones, R2026a Update 4): **421 passed, 0 failed,
0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly: 413 + 8 (`tNgfFitCalibration` multi-arc contract,
reference, precision against a dense joint LS by QR, vectorization, robustness, metamorphic; `tNgfCalibrate`
campaign contract and the three real days). `codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when
re-run after the gate. 59 units (none added: both entry points take arrays).

### N-05 — calibration figures: `app.plotCalibration`, `app.plotAttribution` (issue #62, 2026-10-15)

**Scope, frozen 2026-10-15.** Tier A. CONCEPT 10.1 row N-05 ("calibration figures", the v3 plan's WP-20 part);
timeSeriesPreparation X-11 already draws the run's figures (`report`), so N-05 is the two the calibration adds:
v3 WP-20's `calibration` ("scale/bias with uncertainty bands + residual ACF before/after whitening") and
`attribution` ("per-rev scale vs drivers with regression bands + verdict banner").
- **Afterwards true:** two new units (59), each taking the result objects, never raw arrays, returning figure
  handles, with 'Visible' ("on") and 'ExportTo' (PNG, 200 dpi, `exportgraphics`), as `timeSeriesPreparation.report`.
  `app.plotCalibration(cal, t)`: per axis the scale S_ii(t) with its pointwise band (`scaleAt`, `scaleSigmaAt`)
  and the bias b_i(t) with its band (the basis through the bias covariance); with 'Residual' (fitCalibration's
  `info.residual`) the residual ACF before and after AR whitening (`timeSeriesPreparation.tsmodels.tsAcf`, and
  `tsFitAR(...).resid` for the innovations). `app.plotAttribution(T, A, info)`: the per-revolution scales with
  their sigmas over time, each attributed parameter against its strongest drivers with the fitted partial line,
  and the verdict as a banner. Every drawn object carries a Tag, so a test can read back what was drawn.
- **Out of scope:** an HTML report (X-11's `report` writes one for the run); model-comparison figures.
- **Oracles.** v3 WP-20: "No precision refval (state so)". What a test CAN assert is that the figure draws the
  calibration and nothing else: every curve's data equal the source object's values, bitwise (one authority).
  reference, vectorization and metamorphic are exempt, each with its reason in `testCategoryExemptions.m`.

- **Found while building it, and fixed here.**
  1. **`app.attributeRevolutions` (N-03c, on `main`) refuses fewer than 20 revolutions**: the automatic block
     length of `regressDrivers` needs 20 rows, so ONE day (16 revolutions) raised
     `timeSeriesPreparation:blockLength:tooShort` -- an error its help did not list, beside an ACCURACY line
     claiming a day "supports a few drivers". The help now states the 20-revolution minimum and lists the
     error; a `'BlockLength'` pass-through lets a short arc proceed on the caller's choice; its robustness test
     provokes the refusal and pairs it with the accepted call.
  2. **The whitening is Yule-Walker on the longest gap-free run, applied as a filter to every run.**
     `tsFitAR`'s default exact Kalman method accepts gaps but did not finish 4000 samples inside a minute.
- **Measured:** the calibration figure of a quarter day draws in 1.9 s; the attribution figure of two days in
  0.3 s; with a planted secular scale drift the attribution reads "sensor" for S_xx.

**Delivered 2026-10-15, ☑ on Tier A** (`ciGate` over fresh clones, R2026a Update 4): **413 passed, 0 failed,
0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly: 405 + 8 (`tests/tNgfFigures.m`: contract,
precision, robustness, speed for each unit; reference, vectorization and metamorphic exempt). `codeIssues` on
`+NonGravForces` empty; `auditPackage.py` clean when re-run after the gate. 59 units. Branched from release 0.15.0
(#87), so the PR's diff narrows to N-05 once #87 merges.

### N-03c — per-revolution fit and its attribution on drivers (issue #60, third of three PRs, 2026-10-15)

**Scope, frozen 2026-10-15.** Tier A. CONCEPT 6.2 ("Attribution: per-revolution parameters are regressed on
per-effect a priori amplitudes, space weather, beta' and u; a VIF above 10 gives the verdict unresolved") and 7
(drivers supplied to `timeSeriesPreparation.regressDrivers`, which owns the regression, the VIF gate and the
block-bootstrap sigmas); the v3 plan's WP-16 `perRev` and WP-17 `driverRegression`.
- **Afterwards true:** two new units (57).
  - `app.perRevolution(obs, ref)`: contiguous windows of 'Period_s' (default 5640 s) from the first usable
    epoch; per window `app.fitCalibration` with a constant bias and a diagonal scale (S_yy held by default, as
    there). A window below 'MinWindowFraction' (0.8) of usable rows, or one the estimator's own gates refuse
    (`timeSeriesPreparation:observability:insufficientData`, `NonGravForces:fitCalibration:numericCollinear`),
    is a row of NaN with the reason recorded -- never a number (v3 WP-16). Other errors propagate. Optional
    'Effects' (the per-effect a priori arcs on the raw grid) add each effect's RMS per window and axis.
  - `app.attributeRevolutions(T)`: per-revolution drivers -- from 'Orbit' beta' (`referenceSystem.
    orbitSunGeometry`) and the eclipse fraction (`env.eclipseFactor`, conical), from 'SpaceWeather' F10.7, its
    81-day mean and daily Ap, the effect amplitudes T carries, the window time, and any caller 'Drivers' --
    then `timeSeriesPreparation.regressDrivers` on each estimated per-revolution scale. Each driver belongs to
    a family (density: f107, f107a, ap, amp_drag; geometry: betaPrime, eclipseFraction, amp_radiation,
    amp_thermal; sensor: time); the verdict per parameter is the family of the drivers whose bootstrap
    |coef/sigma| exceeds 2, "mixed" for several families, and "unresolved" with its reason when
    regressDrivers refuses for collinearity (the CONCEPT's VIF verdict) or nothing is significant.
- **Out of scope:** the 1/rev terms and their refused combination with a per-revolution scale (no 1/rev columns
  exist yet); `decomposeScale` (polar decomposition of a full S, v3 WP-17), which the CONCEPT does not list;
  multi-arc campaigns (N-09).
- **Oracles.** perRevolution: a noise-free pair returns every window's truth to the computed rounding bound;
  all windows equal the per-window fits done one at a time (bitwise). attributeRevolutions: the coefficients
  equal `regressDrivers` called directly on the same drivers (bitwise); v3 WP-17's synthetic cases (a
  per-revolution scale built as 0.9 + 0.05 z(F10.7) + noise gives "density" with the F10.7 coefficient within
  20 %; a pure secular drift gives "sensor"); beta' of a circular orbit whose plane holds the Sun is 0.

- **Measured while building it, and acted on.**
  1. **The verdict threshold is 3 bootstrap sigma, not 2.** On 60 synthetic revolutions regressDrivers' block
     bootstrap chose a block length of 14 (about four independent blocks), and the sigma of a driver with NO
     effect came out 0.54 of the OLS sigma: a null driver read |z| = 2.7, and both v3 cases came out "mixed" at
     2. At 3 the verdicts were right in 12 of 12 seeded cases (6 density, 6 sensor).
  2. **Daily indices over a few days make bootstrap replicates rank deficient** inside regressDrivers
     (`MATLAB:rankDeficientMatrix`, 40 revolutions = 2.6 days): f107, f107a and ap are step functions there.
     The point estimate is unaffected; the help says to attribute space weather over weeks (N-09), and the one
     test that provokes it silences that identifier only.

**Delivered 2026-10-15, ☑ on Tier A** (`ciGate` over fresh clones, R2026a Update 4): **405 passed, 0 failed,
0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly: 391 + 14 (`tests/tNgfRevolutions.m`, seven
categories for each of the two units). `codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when re-run
after the gate (its first run caught commas in a catalogue h1 and three size cells written `Rx*`). 57 units.
**N-03 is complete; its release is O-49.**

### N-03b — model comparison by blocked cross-validation, `app.compareCalibrations` (issue #60, second of three PRs, 2026-10-15)

**Scope, frozen 2026-10-15.** Tier A. The v3 plan's WP-14 `compare`: "fit each candidate spec, rank by
blocked-CV score; report AIC/BIC only as reference columns flagged IC_unreliable"; CONCEPT 3.2: no
estimator in timeSeriesPreparation, so `timeSeriesPreparation.crossval` partitions and the application fits.
- **Afterwards true:** new unit `app.compareCalibrations(obs, ref, specs)` (55 units). `specs` is a cell of
  `fitCalibration` name-value cells; ONE `crossval` call scores every spec per fold and per axis (the loss
  function fits all of them on the training rows and returns their test-row mean squared prediction
  errors), so every spec sees the same folds. Bootstrap is off throughout (ranking needs no uncertainty).
  The buffer defaults to the larger of crossval's automatic block length and the CRN window in rows, so no
  test row's filter window shares a raw sample with a training row's. Result: a table, one row per spec,
  with per-axis `cvMse`, `cvMseSe`, `rank`, the parameter count, and `aic`/`bic` from the full-arc fit as
  reference columns whose descriptions say `IC_unreliable` (the residual is CRN-filtered and correlated, so
  the effective sample size is not n); `UserData.best` names the winner per axis.
- **`app.fitCalibration` gains** `'Span'` (int64 [first last] the bases must cover; default the usable epochs,
  as before, so every existing call is bitwise unchanged) and `info.prediction` (the fitted observation at
  every row whose design is finite, fixed scale entries included). Without `'Span'` an edge fold's B-splines
  would be ZERO over the test block, not extrapolated.
- **Out of scope:** N-03c; choosing the number of folds; a CV-selected penalty weight.
- **Oracles.** The composition (one fold re-derived by fitting `fitCalibration` on its training rows by hand,
  bitwise); batched specs equal specs compared one at a time (vectorization); the selection rate over seeded
  realisations against the planted bias order, with AIC's rate recorded beside it.

- **Found while building it, and fixed here.**
  1. **An N-03a defect on `main`:** `fitCalibration` with `'BiasBasis' "constant"` and S_yy held failed with a
     dimension error. Axis y then has ONE parameter, lsqvce returns a scalar, and a scalar indexed by an empty
     index is 1x0 where a column is 0x1, so `J * x(sIdx)` broke; the pre-N-03a `x(sIdx).'` had tolerated it.
     No test fitted a constant bias with the default fixed cross-track scale, so #83's gate could not see it.
     Fixed by forcing the column shape, with a regression assertion in `fitCalibrationRobustness`.
  2. **The minimum-data gate is scaled on a fold.** With 2 folds a training set always keeps under half the rows,
     and `fitCalibration`'s 0.5 weight gate read the held-out block as lost data and refused. On a fold the
     fraction is the candidate's own times the share of finite rows the fold keeps; the full-arc fit is
     unchanged.
  3. **A local basis whose knots do not exceed the held-out window is refused** (`inputFolds`). A B-spline with
     1800 s knots on a 0.25-day arc, 5 folds, has functions whose in-arc support lies wholly inside an edge
     block; lsqvce rightly called that fold singular. Refusing up front beats failing mid-fold, and beats
     skipping the folds a candidate cannot predict, which would score it better for it.
- **Measured.** On the default pair (bias drifting on x and z, not on y), blocked CV picked the planted order on
  9 of 9 axes over three seeds; AIC and BIC picked the 5400 s B-spline on 9 of 9, whose CV error was four orders
  of magnitude worse (it cannot extrapolate over an edge block).

**Delivered 2026-10-15, ☑ on Tier A** (`ciGate` over fresh clones, R2026a Update 4): **390 passed, 0 failed,
0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly: 383 + 7 (`tests/tNgfCompareCalibrations.m`, all
seven categories); `fitCalibration`'s additions extend its existing contract, robustness and metamorphic tests.
`codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when re-run after the gate. 55 units.

### N-03a — time-varying scale in `app.fitCalibration` and `app.Calibration` (issue #60, first of three PRs, 2026-10-15)

**N-03 is three PRs under #60, as N-02 was** (CONCEPT 10.1 row N-03: "time-varying scale, model comparison,
per-revolution, attribution", the v3 plan's WP-14, WP-16 part, WP-17 part): **a** the time-varying scale
(WP-14 `timeVaryingFit`); **b** model comparison (WP-14 `compare`, ranked by blocked CV through
`timeSeriesPreparation.crossval`); **c** the per-revolution fit (WP-16 `perRev`) and its attribution on drivers
through `timeSeriesPreparation.regressDrivers` with the VIF verdict (WP-17; `orbitSunGeometry` from R-01 supplies
beta' and u). The release follows c.

**Scope of N-03a, frozen 2026-10-15.** Tier A. CONCEPT 6.2: "time-varying s(t) only with its pointwise sigma band
and an isDriftSignificant flag".
- **Afterwards true:** `app.fitCalibration` takes `'ScaleBasis'` ("constant" default | "linear" | "hat" |
  "pspline") and `'ScaleKnotSpacing_s'` (43200). Every estimated S_ij becomes S_ij(t) = B_s(t) g_ij on the basis
  `Calibration.basis` evaluates, with columns H(B_s,q .* a_j): the scale multiplies the raw signal before the
  filter. A knotted scale finer than a knotted bias is refused (`inputScaleBasis`). `app.Calibration` carries
  `scaleCurve` (.basis .t0 .knots .tRef .coef 9xQ .cov 1x9), `scale` = S(tRef) checked against the curve at
  construction (one authority), `scaleAt`, `scaleSigmaAt` (the pointwise band), `drift` (Wald test of constancy,
  chi-square, p from base MATLAB's `gammainc`) and `isDriftSignificant`; `apply` solves per epoch as pages;
  `invert` refuses a curve (`unsupportedTimeVarying`); schema "2", and schema "1" files still load.
  `app.syntheticPair` takes `'ScaleDrift'` (3x3 per day), added only when non-zero so every existing pair is
  bitwise unchanged. `app.calibrate` reaches all of it through `'Fit'`, unchanged.
- **Out of scope:** N-03b and N-03c; a smoothness PENALTY on the scale (see the decision below); a time-varying
  K2.
- **Oracles.** Noise-free: S(t) a with a linear S(t) lies exactly in the span of H([1, t] .* a), so the curve is
  recovered to the computed rounding bound the N-02b tests use. Independent: `apply` against a per-epoch loop
  of 3x3 solves with S(t) written by hand; the Wald statistic of the linear basis against slope^2 / var(slope)
  and its p-value against erfc(sqrt(W/2)). Metamorphic: "constant" is bitwise the default path; a hat with one
  span and "linear" are two parameterisations of one column space.
- **Measured while building it, and acted on.** On a one-day pair with a planted 1 %/day drift on S_xx and none
  on S_zz, the FORMAL white-noise covariance (NBoot 0) flagged the undrifted S_zz at p = 1e-18; the bootstrap
  covariance (100 replicates) gave S_zz p = 0.72 and S_xx a slope of 0.01008 +- 0.00074 per day, the bootstrap
  sigma 26x the formal one. So `drift` WARNS (`Calibration:formalDrift`) when its covariance is the formal white
  one: the CRN-filtered residual violates that law by an order of magnitude, and a test on it is not a test.

**Delivered 2026-10-15, ☑ on Tier A** (`matlab -batch` over `ciGate('DependencyRoot', <fresh clones at origin/main>)`, R2026a Update 4): **383 passed, 0 failed, 0 incomplete, 0 shortfalls, GATE TRUE**, the prediction exactly; `codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when re-run after the gate. CI pending on the PR. Prediction for the gate: 375 (after O-47) + 8 = **383**:
`tNgfFitCalibration` +5 (time-varying contract, reference, precision, robustness, metamorphic),
`tNgfCalibration` +3 (time-varying contract, precision, metamorphic); `syntheticPair`'s `ScaleDrift` assertions
extend its existing contract and metamorphic tests. 54 units (none added).

### O-47 — the `qc` forwarders deleted; release 0.14.0 (2026-10-15)

**Scope, frozen 2026-10-15 after PR #81 (release 0.13.0, the one release that carried the forwarders) merged.**
Tier A. CONCEPT 10.1 N-01: "qc forwards for one release, then is deleted."
- **Afterwards true:** `+NonGravForces/+qc/` is gone (six units and its `Contents.m`), with `tests/tNgfQc.m`,
  `tests/tNgfQcScreen.m`, the six catalogue rows, `NgfTestCase/hushQcDeprecation`, `eventTable`'s four
  category exemptions, and the L7 line of every layer map (`CLAUDE.md`, the user guide, the theory body,
  `buildDocs.py`'s `LAYERS`). L7 is left VACANT rather than renumbering `app`, so every "L8 app" written
  anywhere stays true. The hampel and PELT oracles (`tools/pyxcheck/hampel_*`, `pelt_*`) went too: their only
  readers were the deleted suites, and timeSeriesPreparation carries its own copies under
  `tests/tools/pyxcheck/` (read 2026-10-15). Every remaining caller in help text -- `sensor.Electrostatic`,
  `sensor.Sensor`, `io.readTHR1B`, `provenance` -- names `timeSeriesPreparation.qc.*` (or, for provenance's
  example, `app.forceBudget`, since timeSeriesPreparation stamps its own); `tNgfIo/readTHR1BContract` compares
  the THR1B schema with `timeSeriesPreparation.qc.eventTable()`. 60 -> 54 units.
- **Version: 0.14.0.** Removing six public units breaks a caller, which is MAJOR from 1.0.0 on; ci-actions
  `VERSIONING.md` gives `0.y.z` "no stability promise of any kind", so the removal is the next MINOR, as
  0.13.0's replacement of `app.calibrate`'s signature was. `tools/benchmark.m` re-run on Tier A for it
  (audit check 13).

**Delivered 2026-10-15, ☑ on Tier A** (`matlab -batch` over `ciGate('DependencyRoot', <fresh clones of
referenceSystem, gracefoReading, odeSuite, timeSeriesPreparation, lsqvce at origin/main>)`, R2026a Update 4):
**375 passed, 0 failed, 0 incomplete, 0 shortfalls, GATE TRUE** (`tests/last_run.txt`). The first run was 374/375: `tNgfCore/toolboxVersionContract` asserted at least 8 `Contents.m` files, a literal that was the old package count; it now derives the count from the `+pkg` folders (root + one each), which is stricter, and the second run was green. `tools/benchmark.m` was re-run for 0.14.0 (median of 5) on a host running about twenty other MATLAB processes: the compiled `kernel.singleEpoch`, which this change cannot reach, read 0.497 ms against 0.13.0's 0.457, and the MATLAB rows 5-30 % higher, so the shift is the host, not the deletion. The user guide's Cost section and `CLAUDE.md` quote the new table. Prediction for the gate was 414 - 39 = **375**: the 39 is MATLAB's own
`TestSuite.fromFile` count of the two deleted classes (tNgfQc 11, tNgfQcScreen 28), taken before deletion.
`codeIssues` on `+NonGravForces` empty; `auditPackage.py` clean when re-run after the gate.

### N-02c — `app.calibrate` rewritten over `timeSeriesPreparation.run`, with N-06's records (issue #59 third PR, #63; 2026-10-15)

**Scope, frozen 2026-10-15.** Tier B. CONCEPT 4 (the entry `calibrate(acc, M, 'Orbit', orb, 'Attitude', att)`
| `calibrate(acc, refArc)` and its result fields), 6.1 (stages; a jump is a Heaviside column, never a
mask), 6.2 (a priori per effect from the Model; attribution with a VIF gate); N-06 (#63): A.2 admits
the app layer, the dependency table carries timeSeriesPreparation and lsqvce, B.4 records the re-entry.
The 0.12.0 matrix interface `calibrate(yMeas, aModel, epoch)` and its own QR are REPLACED outright
(owner ruling 2026-10-14).
- **Afterwards true:** `app.calibrate(acc, reference)` takes an observation Arc and a Model (with
  'Orbit', optional 'Attitude') or a raw-grid reference Arc. A Model is evaluated along the orbit per
  effect and interpolated by pchip onto the CRN kernel's raw grid; `timeSeriesPreparation.run` runs
  mask, crn, align and qc; qc's "jump" events (UTC) become GPS-ns bias steps; `app.fitCalibration`
  estimates; a per-arc attribution regresses the residual on the filtered per-effect series through
  `lsqvce.fit` with the VIF gate. `res` is run's struct plus `calibration` and `attribution`
  (`diagnostics.fit`, `diagnostics.rawReference`). `app.Calibration` carries jump steps
  (`bias.jumps`); `fitCalibration` takes 'Jumps' and 'Robust' and drops, with a warning, a jump whose
  step is constant over the fitted rows. `syntheticPair`'s observation is labelled ACT1B. referenceSystem's
  minimum rises to 2.3.
- **Out of scope:** the per-revolution attribution on drivers (N-03), multi-arc (N-09), the release
  itself (O-48: the version bump needs `tools/benchmark.m` re-run on Tier A, which Tier B cannot).
- **Decisions for the owner (Part I, 2026-10-15):** the A.2 rule's new wording; the Model evaluated
  at the orbit's epochs and interpolated, not at 10 Hz; jumps shared by the three axes.
- **Measured on the first CI run (033ac76, cancelled by the 45-minute job timeout), and changed:** with
  qc's jump events fitted by default, `tNgfCalibrate/calibratePrecision` ran 16.8 minutes and put the
  scale 22.4 bootstrap sigma off the truth, and the calibration demo ran into the timeout. qc screens
  the UNCALIBRATED residual, where a scale error and a bias drift look like level changes to PELT, and
  every spurious jump became a parameter. `'Jumps'` is therefore OPT-IN (default false); the epochs are
  still converted and reported in `diagnostics.jumps`. The same run measured `calibrateSpeed`'s 2x-arc
  ratio at 4.26 against a bound of 4 whose premise ("every stage linear") was false -- qc's change-point
  search is quadratic in the worst case; the test now asserts calibrate's own overhead against its parts.

**Delivered 2026-10-15; ☑ on the CI runner** (run 38704980955 on the merge commit of #80: 414 passed, 0 failed, 0 incomplete, 0 shortfalls, GATE TRUE -- the prediction below, exactly). Prediction for the gate: 414 (after N-02b) - 7 (the old calibrate
tests in `tNgfApp`) + 7 (`tests/tNgfCalibrate.m`: contract, reference, the real-day integration,
precision, robustness, metamorphic, speed; vectorization exempt) = **414**, 0 failed, 0 incomplete,
0 shortfalls, 60 units. The reference oracle is the composition (bitwise equal to fitCalibration on
run's arc); the integration runs ACT1B of 2020-01-01 against the fast Model along GNV1B with SCA1B.

### N-02b — `app.fitCalibration`, the estimator on `lsqvce.fit` (issue #59, second of three PRs, 2026-10-15)

**Scope, frozen 2026-10-15.** Tier B. CONCEPT 6.2 and 6.6 (item 14): per axis, y_i = H*b_i +
sum_j S_ij H*a_j [+ K2_i H*a_i^2] + e_i, every deterministic column filtered by the same CRN kernel as
the observation, estimated by `lsqvce.fit`; the three axes are three problems. `app.calibrate` is
untouched until N-02c replaces it; the new estimator is its own public unit so that N-02c's entry
(arc or `Model` reference, `run`, attribution) can wrap it without re-opening the numerics.
- **Afterwards true:** `app.fitCalibration(obs, ref)` returns an `app.Calibration`. Bias basis constant |
  linear (default) | hat | pspline on knots spanning every raw sample a window reads (warns below the
  kernel span); scale diag (default) | full | mask, S_yy held at its a priori by default; K2 optional.
  Kind A per axis: one system, one day, one arc per contiguous run of usable rows; columns equilibrated
  to unit RMS. Gates from `timeSeriesPreparation.observability` (minimum-data gate, cond, paramCorr,
  VIF); a row with several model columns and a VIF above 'VifMax' (10) is REFUSED. Uncertainty: formal
  (lsqvce Cxx) and a moving-block pairs bootstrap through `timeSeriesPreparation.bootstrap`, refitting
  every replicate with `lsqvce.fit`; the Calibration carries the bootstrap covariance. EIV bound per
  model axis (v3 D4, Fuller 1987). lsqvce re-enters `dependencyReport` (optional, minimum 6.22).
- **Out of scope:** jumps and 1/rev columns (they need `run`'s QC and the orbit, N-02c), the
  per-revolution and time-varying scale (N-03), multi-arc (N-09), N-06's prose records (N-02c).
- **Deviation from the concept text, argued (RATIFIED by the owner 2026-10-15):** the K2 column is H(a_i^2), not (H a_i)^2 as CONCEPT
  6.2 writes it, because the nonlinearity acts on the instrument's raw signal before the filter.
- **Decisions taken here, for the owner to overrule:** see Part I, 2026-10-15.
- **Measured on the first CI run (38690853810), and kept as evidence for the white default:**
  `'Noise', "estimate"` on a 0.05-day synthetic arc ran lsqvce's joint estimator to its 10-pass budget
  without settling (last max|alpha^2 - 1| = 6.4, 6.8 and 6.4 on the three axes): the law estimated
  from a CRN-band-limited residual does not converge. The same run found two defects of this PR's
  own: the cofactor was scaled back with the observation scale (a bound 5e6 too tight), and the
  knots spanned the whole arc (a hat outside every fitted window, refused as singular). Both fixed.

**Delivered 2026-10-15, ☑** on the CI runner after the two fixes (run of f9587d4: **414 passed, 0 failed, 0 incomplete, 0 shortfalls, predicted 414**). Prediction for the gate was: 408 (after N-02a) + 6 = **414**
(`tests/tNgfFitCalibration.m`: contract, reference, precision, robustness, metamorphic, speed;
vectorization exempt with its reason), 0 failed, 0 incomplete, 0 shortfalls, 60 units. The
noise-free reference recovers a full S with off-diagonals, the drifting bias and K2 against a bound
computed per parameter as 2 sqrt(Q_kk n) eRow + 100 eps cond sqrt(Q_kk n) rms y.

### N-02a — the calibration truth and its result object (issue #59, first of three PRs, 2026-10-14)

**Scope, frozen 2026-10-14.** Tier B. The owner split N-02 into three PRs under #59 (decision log):
N-02a the truth generator and `app.Calibration` with save/load; N-02b the estimator on `lsqvce.fit`;
N-02c the arc/model entry, `timeSeriesPreparation.run`, attribution, the integration test, N-06's records
and release 0.13.0, which replaces `app.calibrate`'s signature outright. Until N-02b, `app.calibrate` is
untouched and nothing returns a `Calibration`.
- **Afterwards true:** `app.syntheticPair` builds a 10 Hz model arc and passes S·model + b(t) + K2 +
  coloured noise + 1 s spikes + one gap through the REAL L1B chain (`crnKernel("ACC")`, `crnApply` at the
  whole seconds), returning every planted number; all randomness comes from `'Seed'`. `app.Calibration`
  holds the observation-direction parameters (y = b(t) + S a), a bias CURVE on a named basis (constant,
  linear, hat, cubic B-spline on int64 knots) whose evaluator `Calibration.basis` is the one the N-02b
  estimator will build its design from, `apply` (refuses a K2 term rather than invert it approximately),
  `invert` (TN-02's correction direction, first-order covariance, documented as not the reverse
  regression), `summaryTable`, and a versioned MAT file (`schemaVersion` "1"; a newer schema warns).
- **Out of scope:** any estimation (N-02b); `app.calibrate` (N-02c).
- **Found on the way:** `tools/audit_class.py` and `tools/audit_args.py` split the name of a function
  WITHOUT outputs: `(?:...|\w+)?\s*=?` let `function save(obj, file)` read as a method `e`, so the help
  check of every output-less public method had been run against the wrong name. The output clause is
  now optional as a whole. No existing unit was affected (the audit stayed clean after the fix).

**Delivered 2026-10-14, ☑** on the CI runner (run 38685265825, PR #78: **408 passed, 0 failed, 0 incomplete, 0 shortfalls, predicted 408**). Prediction for the gate was: 394 (after N-04 and N-01) + 14 = **408**
(`tests/tNgfCalibration.m`, seven categories for each of the two units), 0 failed, 0 incomplete,
0 shortfalls, 59 units. Every bound is computed from Nf·eps·Σ|F|·max|raw| or cond(S)·eps; the
B-spline basis is checked against the closed-form uniform cubic, not against Cox-de Boor.

### N-01 — `qc` forwards to timeSeriesPreparation (issue #58, 2026-10-14)

**Scope, frozen 2026-10-14.** Tier B. CONCEPT 10.1 N-01: "qc forwards for one release, then is deleted;
O-1 and O-18 move with it." timeSeriesPreparation's X-05 ported the six units unchanged and X-16 made them
fast under an isequaln gate; the two argument blocks were compared and are identical.
- **Afterwards true:** each `NonGravForces.qc.<unit>` warns `NonGravForces:<unit>:deprecated` and calls
  `timeSeriesPreparation.qc.<unit>` with the same arguments; without the target it raises
  `NonGravForces:<unit>:depMissing`; the target's errors propagate unwrapped (`timeSeriesPreparation:*`).
  `io.readTHR1B` builds its table with the target directly. The pre-move suites `tNgfQc` and
  `tNgfQcScreen` run unchanged through the forwarders, only the identifiers updated. O-1 and O-18 close
  here as moved; O-47 opens for the deletion.
- **Out of scope:** the deletion itself (O-47, after the release that carries the forwarders).

**Delivered 2026-10-14, ☑** on the CI runner (run 38636885801, PR #77 with N-04: **394 passed, 0 failed, 0 incomplete, 0 shortfalls, predicted 394**). Prediction for the gate was: 393 (after N-04) + 1 = **394**
(`tNgfQc/qcForwardingContract`: bit-equality with the target, the warning provoked on all six units,
`depMissing` provoked by taking the package off the path), 0 failed, 0 incomplete, 0 shortfalls.

### N-04 — `sensor.Electrostatic` takes its default kernel from `timeSeriesPreparation.crnKernel` (issue #61, 2026-10-14)

**Scope, frozen 2026-10-14 before any code.** Tier B. Owner ruling 2026-10-14 (asked, three options
put): the DEFAULT kernel is `crnKernel("ACC")`, applied to 10 Hz dynamics; a 1 Hz arc is refused unless
the caller passes `[]` explicitly. Plumbing went first as its own PR (#76: timeSeriesPreparation and lsqvce
as gate siblings).
- **Afterwards true:** `Electrostatic(bias, S, dr)` filters with the toolbox kernel's taps, unchanged;
  `'SampleRate'` defaults to 10 and must equal a CRN kernel's `fs`; `dyn.epoch` (seconds, datetime or a
  referenceSystem `.day/.sec` epoch) is required with a CRN kernel and checked to lie on its grid;
  `'OutputRate'` samples the filtered series (never averages) as ACT1B is formed; `dependencyReport`
  lists timeSeriesPreparation (optional, minimum 1.4); docs, catalogue, user guide and `CLAUDE.md` follow.
- **Out of scope:** using `crnApply` itself (it needs a timeSeriesPreparation `Arc` on int64 ns epochs;
  the sensor convolves the same taps with the same sum), a version bump (planned 0.13.0 with N-02/N-06),
  and the calibration upgrade (N-02).

**Delivered 2026-10-14, ☑ (Tier B; verified with N-01 by the CI `matlab` job, run 38636885801: 394/394, predicted 394).**
- **The oracle is the published number, recomputed:** the DC gain through the sensor is asserted against
  the printed 1 + 1.115e-7 (±5e-11, half the last printed digit), and an independent Python evaluation
  of ATBD App. B written for this check gives 1 + 1.11516e-7 and a gain at f0 of exactly 1. Unit gain at
  f0 and zero phase (the kernel is symmetric) are asserted on a sinusoid, with the bound computed from an
  Nf-term dot product's rounding.
- **Refusals are provoked, each beside a call that must succeed:** no `dyn.epoch`, a 1 Hz grid, a wrong
  length, an unsupported epoch type, a kernel rate other than `'SampleRate'`, text other than "crn", a
  struct without `.F/.fs`, a non-integer output ratio; and `depMissing`, by taking the package off the
  path for one call and restoring it.
- **Seventeen existing sensor constructions now pass `[]`**, because they test unfiltered 1 Hz physics; the
  `app.calibrate` help example does the same with `'SampleRate', 1` (1 Hz L1B is already filtered).
- **Superseded decision:** "The CRN filter is an argument, not an implementation" (2026-09-08) held while
  no kernel had been read from the ATBD; timeSeriesPreparation's `crnKernel` now has, against an
  independent transcription (its ORC-1) and real ACT1A/ACT1B (its X-12). The argument remains; the
  default is no longer empty.
- **Prediction for the gate:** 388 + 5 = **393** (tNgfSensor + 5: `ElectrostaticCrnContract`, `...Reference`,
  `...Precision`, `...Robustness`, `...Metamorphic`; tNgfCore unchanged in count), 0 failed, 0 incomplete,
  0 category shortfalls.
- **Next:** N-01 (#58), `qc` forwards to timeSeriesPreparation's port.

### N-07 — `geom.grace`, the GRACE-A/B macro model (issue #64, opened 2026-09-22)

**Scope, frozen 2026-09-22 before any code** (timeSeriesPreparation GUIDE step 3). Tier B: no MATLAB in
this session; the MATLAB suite runs in CI on the self-hosted runner.
- **Afterwards true:** `NonGravForces.geom.grace("A"|"B")` returns a `MacroModel` whose every number is
  read from R-5 (Bettadpur 2012, GRACE 327-720 Rev 4.6, Sec. I.5.3.3-I.5.3.5, PDF pp. 20-23, supplied by
  the owner on 2026-09-22 because the session's network reaches no host that serves it); tests in the
  seven categories; catalogue row, `Contents.m`, help, generated HTML and LaTeX, `CLAUDE.md` and `README.md`
  in sync.
- **Out of scope:** classic calibration presets, thermal handling and the mission mass (N-08); a
  `cylinder` element for the boom (no kernel consumes one); any number R-5 does not print (inertia,
  CoG -> proof-mass offset). No version bump: the MINOR bump for the added unit belongs to the release
  that carries it (GUIDE step 13).

**Delivered 2026-09-22, ☑ — verified on the self-hosted CI runner, run `35771831961` (commit `214f632`): `GATE TRUE (388 passed, 0 failed, 0 incomplete, 0 audit shortfalls)`, 14 exempt units. Predicted 388, measured 388. Tier A (local gate) not run: no MATLAB in the session; the runner's run is the MATLAB verdict. The run log of the merge head (run `35773012179`, R2026a Update 5, 388/388) is committed as `tests/last_run.txt`; the owner supplied the artefact, which this session's network could not fetch.**
- **R-5 prints the same surface table as D-56935**, digit for digit, down to note 2's 0.65/0.72 and the
  1.07 budget the non-operating value produces. `geom.grace` transcribes it from R-5 independently and
  `tNgfGeom/graceMetamorphicAcrossDocuments` compares the two transcriptions: one fact, two documents, now
  compared.
- **The geometry is R-5's own**: 3123 / 720 / 694 mm (GRACE-FO: 3122.5 / 721 / 695), and the origin is
  Fig 3's C.o.G., 275 mm above the nadir plane rather than at mid-height; the thruster table confirms it
  (A14/A22 at Z = +275 mm, A12/A24 at Z = -444 mm, O11/O21 at X = -1561 mm). The 2 cm non-closure of the
  printed dimensions is the same as GRACE-FO's and is resolved the same way.
- **Checked against numbers the unit does not use** (measured in Python before any MATLAB run): the
  overhang reaches 55.998 mm below the nadir plane (Fig 3: 56), and the side silhouette plus the printed
  boom area is 2.4696 m^2 against Fig 4's 2.4633921 (+0.25 %; the dome, notch and recesses are not plates).
  Polygon areas are within 2.14 % of the printed ones (front and nadir -2.14 %, arrays +0.26 %).
- **Not in R-5, not set**: the inertia tensor and the CoG -> proof-mass offset are NaN, so the
  gravity-gradient torque is refused (`kernel.bodyTorques:inputMacro`) instead of guessed. The boom is
  excluded for the same reason as in `graceFO` (no kernel consumes a cylinder), though R-5's note 4 does fix
  its axis along Z.
- **Found while reading R-5**: classic MAS1B is spacecraft mass, not tank mass (O-46).
- **Prediction for the gate:** 381 + 7 = **388** tests (tNgfGeom: 10 + 7 = 17), 0 failed, 0 incomplete,
  0 category shortfalls (`grace` holds exemptions for vectorization and speed only, as `graceFO`).
- **Read from the run log, 2026-09-22:** every `grace` bound held, and the measurements match the Python
  predictions (side silhouette 2.5328e-3 vs 2.53e-3; overhang 2.2 um from 56 mm). One bound was a knife
  edge: `grace zenith plane vs thrusters` measured **1.0000e-03 against 1.0000e-03**, passing on its 1e-12
  only, because the gap is 1 mm by construction. It now asserts the thruster lies inside the face and bounds
  the gap at 2 mm, the 1.5 mm that millimetre rounding of three printed values allows plus margin.
- **Next session:** N-08 (#65) when N-02 and X-13 are in.

## Part H — Obligations and open items

| # | Item | Status | Owner | Re-entry |
|---|---|---|---|---|
| O-51 | **Classic GRACE Level-1B for 2008-01-01 is PARKED (owner, 2026-10-26): GFZ's ISDC is down.** It is the only anonymous source; `isdc-data.gfz.de` answered neither GitHub-hosted runners (runs 40358650063, 40361879447) nor the self-hosted runner (40365636347, 40375328241: WinError 10060), and `isdcftp.gfz-potsdam.de` no longer resolves. No archive name for the classic day has been read, so nothing is pinned. N-08 (#65) needs the A/B ACC1B, SCA1B, GNV1B, THR1B and MAS1B excerpts | open | owner | when ISDC answers: add the day to `setupNonGravForces('FetchPool')` (the one network authority) and cut it in `tools/cutFixtures.py` with the `data/samples` recipe |
| O-52 | **DTM2020 has no native evaluator (2026-10-27).** `app.DynamicsSource` refuses `drag.density = "dtm2020"` at construction (`unsupportedOption`), as it does `ceresClimatology` (O-44): `mex/ngfRhs.c` links the density models it evaluates and CNES's subroutines were not added to it, and an unknown name would otherwise pack as `exponential`. Also open: the RESEARCH DTM2020 (F30, Hp60), which needs the hourly Hp60/ap60 index (GFZ; not in hand) and whose shipped `LECDTM` cannot be called without editing CNES's file, which the licence forbids | open | session | when an integrated arc needs DTM2020: link the subroutine file into the rhs target and pass the four drivers through the tables; the research model when Hp60 is in the pool and CNES ships a callable reader |
| O-53 | **CLOSED 2026-11-03 (owner: change to the new API).** `app.fitCalibration` (joint and single-arc solves, bootstrap refit) and `app.calibrate`'s attribution declare a `lsqvce.TaskDefinition` (`fromMatrices` for kind B; one system "acc", commons `c1..cp`, day-scoped locals `b1..bP` placed by name, stream "acc" with `lsqvce.Noise.white(1)` or the user's law, `Estimate="law"` for "estimate"); `'Noise'` takes an `lsqvce.Noise` and refuses the v6 struct (`fitCalibration:inputNoise`); `dependencyReport`'s minimum is 7.0 and `fitCalibration:depMissing` fires without `lsqvce.TaskDefinition`; CI back on lsqvce `main` -- WHAT BROKE: **lsqvce 7.0.0 broke `app.fitCalibration` (2026-11-03).** lsqvce's `main` released 7.0.0 (then 7.0.1), whose `lsqvce.fit` takes an `lsqvce.TaskDefinition` and refuses the v6 data struct `fitCalibration` builds (`lsqvce:fit:structDeclaration`); CI 41691731161 on PR #132 (documentation only) had 60 calibration tests error, 451/511 passed. `ci.yml` and `release-checks.yml` now check lsqvce out at `v6.35.1`; `dependencyReport` checks only the minimum 6.35, so a v7 on a user's path still fails at the first fit. | closed | Claude | — |
| O-50 | **CLOSED 2026-10-15: setup downloads time out after 60 s, not MATLAB's 5 s, and retry once.** CI run 38786016017 (PR #90) lost NRLMSIS 2.1 to "timed out after 5.000 seconds" from map.nrl.navy.mil -- `setupNonGravForces` set no `weboptions` timeout -- and run 38733381139 (PR #85) a NAIF transfer cut mid-file. All six `websave` calls now go through `localDownload`: `weboptions('Timeout', 60)` (a slow server's first byte, not a transfer limit), one retry after 5 s, and on two failures the error names both causes, which each caller REPORTS as before. Content checks (SHA-256 pins, first lines, GROOPS magic) are unchanged and run after every download, retried or not. Fault-injected on a verbatim copy of the helper: a refused connection reported both attempts after 9.6 s; the IGRF file (42 kB) arrived in 1.3 s. Tier A gate 421/421, predicted 421 (no test changed) | closed | Claude | — |
| O-1 | **MOVED 2026-10-14 with `qc` (N-01) and closed here.** The unit now forwards to timeSeriesPreparation, which ported it unchanged (X-05) and made it fast under an isequaln gate (X-16); the obligation lives there. Original entry: `qc.hampel` loop is 382× slower than `movmedian` (0.64 s/day, acceptable); first candidate for the density/MEX session or a two-pass movmedian formulation with exact window-MAD | closed | LP-5 | — |
| O-2 | **Closed 2026-09-07 (LP-5):** `MacroModel` packs unit normals (`N = N ./ |N|`) and keeps the printed value in `panels(k).normal`; no kernel renormalises. Measured before the change: thermalRecoil per-plate magnitude off by 7.3e-7 with the printed normals, 1e-16 after. Handbook normals are printed to six digits (norm 1 − 7.3e-7, measured 2026-09-07; the earlier rows declared 5e-7 for the squares); `MacroModel` still refuses beyond 1e-5 | closed | LP-5 | — |
| O-3 | CERES SYN1deg variable names (B7) | open | owner supplies a file | before WP-27 **JUSTIFIED 2026-09-17 by [V], which this row previously lacked.** Vielberg & Kusche build their extended Earth-radiation model on HOURLY CERES SYN1deg and report that using the Earth's outgoing fluxes directly improves the force model, RMS fits down 5 % and 27 % for their two calibrations; their Table 2 grades the source in the order this project has implemented it (Knocke, EBAF albedo/emission, EBAF flux, SYN1deg flux). So the diurnal product is worth the work and someone has measured how much. Note the 1 degree GROOPS climatology delivered in WP-E2 is NOT this: a monthly mean has no diurnal cycle whatever its file is named |
| O-4 | **CLOSED 2026-09-16 (LP-C WP-C9), by removal and not by an adapter.** odeSuite delivered `gjToolbox.DynamicsSource` and `csrc/gj_native.h` (ABI 1); `NonGravForces.app.DynamicsSource` subclasses it, and `Model.rhs` is deleted — a handle matching no signature the integrator takes was worse than no handle at all, and `tests/tNgfModel.m` now asserts the method does not come back. `Model` was not changed otherwise, as this row required. Issue #35 closes with it. Original entry (2026-09-08, narrowed by reading the repository instead of waiting for it): odeSuite was at L0 with no dynamics interface, its validator documented `@(t, r, v)` with numeric t, and `Model.rhs()` returned `@(epoch, r, v, q)`, which would have passed the validator's degraded branch and then been called with the wrong signature | closed | LP-C | — |
| O-5 | `referenceSystem` MEX reported "present but not loadable" by its own setup on 2026-08-27; **2026-09-06: setup reports `mex current` (gateway built 2026-08-05, SOFA 20231011, OpenMP 32 threads), so it loads on the host**; open only for the CI runner. **The LP-6b re-entry is void: that MEX is not linked (B.1), so nothing in this toolbox depends on it loading.** Re-entered only if a C caller ever needs `rsApiRotation` | open | owner | if a C integrator is written |
| O-6 | **PARKED BY THE OWNER 2026-09-18: "O-6 remains open until further notice." It stays open as a recorded gap, and no session should re-raise it, schedule it, or treat it as blocking LP-R until the owner says otherwise. Do not ask again.** The analysis below stands and is what a future session needs if the ruling is lifted. **Exact TRP oracle. Narrowed 2026-09-17 (WP-E5), still open.** [S] Table 3 gives RMS per axis for GRACE-A 2008 and is gated. [H24] was expected to supply the GRACE-FO equivalent and does NOT: its 0.148 nm/s^2 is its own thermal model against its own thermistors on the solar arrays at a fitted (C, k), not thermal recoil against truth, and its Table 3 is a ladder of cross-track RESIDUALS by modelling choice rather than a TRP acceleration. Two checks were taken from it and are gated: case (a) omits thermal emission and leaves 4.22 nm/s^2, which bounds the magnitude of the omitted term for THIS satellite (measured here 5.263 nm/s^2 cross-track, bracketed by 4.22 GRACE-FO and 6.4 GRACE-A), and the (b)->(c) step of 1.34 nm/s^2 says what thermal inertia is worth (measured 0.950 nm/s^2 as a model difference). **What is still missing is a published per-axis TRP acceleration for GRACE-FO** | open | owner | PARKED: only when the owner lifts the parking of 2026-09-18 |
| O-7 | **CLOSED 2026-09-16 (LP-D WP-D4): the kernel is fetched at setup and the path is EVALUATED, not merely refused.** `setupNonGravForces` gains `'FetchKernels'` and pulls `de440s.bsp` (32.7 MB) through `referenceSystem.downloadKernels` - that function calls itself the single authority for where these files come from, so this toolbox repeats no URL and pins only the SHA-256 of what arrives (`addede59...164703`). It is a SEPARATE switch from `FetchData` on purpose: `ciGate` passes `FetchData "never"` so the gate does not test the network for OBSERVATIONS, and an ephemeris is model data like the MSIS and HWM14 sources. **The finding that made this worth doing was the process one, not the numbers**: `rs.de440` is the `precise` preset's OWN DEFAULT and no test had ever evaluated it - the only test naming it asserted its refusal, and the native-parity test substituted the analytic Sun. Now measured: DE440 against referenceSystem's analytic Sun is **0.0031 arcsec** in direction, 0.3 km in range, **4.0e-9 relative in solar radiation pressure** - three orders below the 1e-3 uncertainty of the solar constant itself, so the analytic Sun costs nothing an accelerometer could see, and `precise` now runs as shipped through the C evaluator too | closed | LP-D | — |
| O-8 | **CLOSED 2026-09-09: `readMAG1B` is delivered.** The owner supplied the GFZ ISDC path, the 2020-01-01 archive was fetched, a 600-record excerpt was cut by the same recipe as the other fixtures, and the reader is written against a measured specification. Original entry: MAG1B is reachable and always was. The owner supplied the GFZ ISDC path after my probe went to the wrong host; the 2020-01-01 archive `setupNonGravForces` now fetches holds **65 files**, among them `MAG1B_2020-01-01_C_04.txt` — and also ACT1B, MAS1B, THR1B, TNK1B, IMU1B, AHK1B, KBR1B, GPS1B. There is an oracle for `readMAG1B` whenever someone writes it, and the same archive removes the data excuse from anything else that wanted a real L1B product. What remains is the work, not the obligation | closed | — | — |
| O-9 | Dst and F30 are in none of the three space-weather sources in hand (CelesTrak SW-All; SET SOLFSMY, DTCFILE); DTM2020 (F30, Kp) and any Dst-driven model need another source. **2026-10-26 (PR #107): the data is in hand, no reader is.** Kyoto FINAL Dst for 2007-12, 2008-01, 2019-12, 2020-01 and CLS's daily F30 (gaps interpolated and flagged; Nobeyama's own file has no 1 GHz value on 2020-01-01..03) are committed under `data/spaceweather` with their source SHA-256s, fetched into DATAPOOL by `setupNonGravForces('FetchPool')` and cut by `tools/cutFixtures.py`, both run on the CI runner by `.github/workflows/fetch-data.yml`. Classic GRACE Level-1B (N-08) is parked: GFZ's ISDC is down (owner, 2026-10-26) **CLOSED 2026-10-26 by PR #108**: `io.spaceWeather` reads both (families "dst" and "f30"; 9999 and -1 become NaN; F30 carries the previous day and the centred 81-day mean DTM2020 takes), six tests in `tests/tNgfIo.m`, green in CI 40411254932. What remains is a CONSUMER (DTM2020, LP-3b), which is not this row | closed | LP-3b | — |
| O-10 | **CLOSED 2026-09-07: fixed upstream and verified here.** referenceSystem 2.1.216 trims pre-1972 rows and warns instead of refusing; the whole shipped C04 series now works through `readGNV1B` (measured: 3652 of 23562 rows ignored). Reported as astrobayer/referenceSystem#15 (reproduced: 3652 of 23562 rows precede 1972; a subset from MJD 41317 works). `referenceSystem.interpolateEop` / `rotationItrsToGcrs` refuse the shipped whole `eopc04_20.1962-now` (pre-1972 rows; identifier `leapSecondTable:beforeUtcEpoch`): reported to the referenceSystem owner; NonGravForces documents the constraint (B.1, `readGNV1B` help) and lets the foreign error propagate | closed | owner | — |
| O-13 | `srp.eclipse = "solaarsCF"`: **closed 2026-09-07.** The owner allowed downloads; the open-access thesis (VTechWorks 10919/73537) was fetched to `doc/references/`, Table 3.2 transcribed by three agreeing text extractions into `tools/pyxcheck/solaars_cf_coefficients.csv`; checked by the exact limits (0 and 1), monotonicity, C¹ smoothness and the sigmoid centres of Fig. 3.28. No published numeric ν value exists; the thesis states < 3 % residual against full SOLAARS | closed | LP-3 | — |
| O-14 | MEX build on the CI runner: `buildMex` locates the MinGW-w64 support package itself when `MW_MINGW64_LOC` is unset (measured 2026-09-07: the host had the package at `%ProgramData%\MATLAB\SupportPackages\R2026a\3P.instrset\mingw_w64.instrset` and MATLAB did not know); `ciGate` builds before the suite; a runner without the package turns the gate red by design. **Measured 2026-09-07: the first CI run of LP-3 (PR #11, run 34065197914) built the MEX and went green (11m35s)**: the runner has the package; kept as the rule, no longer open | closed | LP-3 | — |
| O-15 | **CLOSED 2026-09-16 with O-7**: the kernel is available to the suite on every tier, `tNgfEnv/radiationSourcesDe440Reference` evaluates the furnished-kernel path against DE440 itself, and the `precise` parity configuration in `tNgfIntegration` no longer substitutes the analytic Sun. The path is ☑. **A second dependency came with it and is now declared rather than assumed**: NAIF MICE. It is not a repository anyone can check out - compiled MEX, managed centrally at most sites, which is why `referenceSystem` detects it and never adds it to the path - and it was reachable here only because it sits on this machine's saved pathdef. That is the same shape as the two defects of 2026-09-16 (an environment variable and a hand-built MEX the runner lacked), so `dependencyReport` lists it and `ciGate` REQUIRES it by name | closed | LP-D | — |
| O-12 | **CLOSED 2026-09-08.** `gracefoReading` 1.0.1 (its PR #7, issue #6) makes `readL1bMulti` return the same clean empties `readL1b` returns for an all-excluding `TimeRange`, reporting it with `gracefoReading:readL1bMulti:noRecordsInRange`. `io/private/l1bLoad` now routes **every** read through `readL1bMulti`, one file or many; the split and the two hand-fabricated plural aux fields are gone, and `dataEmpty` is reachable for any number of files (`tNgfIo/readACC1BMultiEmptyRobustness`, the path that could not be tested before). The declared minimum moved to 1.0.1 | closed | — | — |
| O-16 | **CLOSED 2026-09-17 (WP-E3).** SESAM (Pilinski et al. 2013) is implemented as `drag.accommodation = "sesam"` on both the MATLAB and the C path, against an oracle that integrates eq. (26a) numerically rather than repeating eq. (26b). What is NOT implemented, deliberately: that paper's DSMC lookup tables and tumbling-sphere shape factor, which map alpha to a drag coefficient -- this toolbox has its own panel geometry and Sentman kernel and needs only alpha. The sticking law was thought to be missing until the paper was read on 2026-09-17; it is eqs (22) and (26b) | closed | LP-E | — |
| O-17 | CLL (Cercignani-Lampis-Lord) and Schamberg gas-surface models: closed flat-plate forms (Walker, Mehta & Koller 2014; Schamberg 1959) and their oracles are not in hand; `gsiCoefficients` refuses them with `unsupportedModel` | open | owner | when a below-300-km use case appears |
| O-18 | **MOVED 2026-10-14 with `qc` (N-01) and closed here.** The unit now forwards to timeSeriesPreparation, which ported it unchanged (X-05) and made it fast under an isequaln gate (X-16); the obligation lives there. Original entry: `qc.detectJumps` PELT loop costs ~0.12 ms per sample in MATLAB (measured 2026-09-07: 2.4 s for 20000 samples; inlining the segment cost changed nothing, the per-iteration interpreter overhead dominates), i.e. ~10 s per day at 1 Hz and the whole cost of `qc.screen`; candidate for the vectorised candidate bookkeeping or a C layer — note that the LP-6b C layer as delivered ports the FORCE kernel only and contains no qc, so it does nothing for this | closed | LP-7+ | — |
| O-20 | `gsiCoefficients` evaluates the leeward moments by the direct erfc difference, which cancels: measured by the mirror's 80-digit reference, ~9e-12 relative at s = 15, theta = 112 deg (x = -5.7) against 2.9e-15 for the mirror's continued-fraction form. The affected coefficients are 1e-14 of the head-on value and carry no force, so nothing is changed in the frozen kernel; the mirror records the better-conditioned form (erfc continued fraction descended two levels, a Veltkamp/Dekker split of u^2, and exp(-u^2 + log f) instead of a product of a subnormal) if a future level needs it | open | LP-6b | if a leeward-dominated regime ever appears |
| O-21 | `radiationPanelForce` classifies a plate as lit by the sign of a naive dot product. The mirror's reconciler found by exact rational arithmetic that at an exact terminator the true n.d can be +4e-17 while the double dot product returns +0.0, so the plate is called dark. The force at stake is 1.8e-22 N. Documented, not fixed | open | LP-6 | if a terminator-crossing study needs it |
| O-28 | **CLOSED 2026-09-17, done in LP-C WP-C1 (finding F2).** `geom.thermalState` takes `'Exposed'` as **NxKxM** -- per plate, per epoch AND per source -- not the single fraction per panel and epoch this row described. The defect it warned about was the one measured on 120 real epochs: the nadir panel read 0 K at 11 of them because the Sun's shadow column zeroed its Earth heating. Found stale by the same sweep as O-23. | closed | — | — |
| O-29 | **CLOSED 2026-09-18 -- and it had already been DISCHARGED on 2026-09-16, by WP-C3, which this row commissioned.** G.2 records the delivery in plain words, in figures themselves now superseded: "WP-C3 cost: 1.21 ms per epoch in the standard preset, 0.30 ms in `fast` ... against O-29's 88 ms". The cube map, the on-demand provenance and the single environment evaluation did exactly what this row asked for, the claim at G.1 ("at least ten times cheaper") was met, **and nobody closed the row.** So it went on saying 88 ms for two more days, and six other places went on quoting it -- `CLAUDE.md` twice, `demos/demoForceBudget.m`, `tests/tNgfApp.m`, `tests/tNgfMex.m` and, on 2026-09-17, a user guide chapter written fresh from it and shipped in the PDF -- while `tests/tNgfModel.m` had held the right value since 2026-09-15. **Re-measured 2026-09-18 by `tools/benchmark.m`, whose output is committed at `tools/benchmark_reference.csv`** (0.12.0, R2026a, median of 5 after five warm-up calls): standard preset **15.9 ms at K=1, 2.99 at K=10, 0.798 at K=100**; fast preset 4.00, 0.695, 0.076; `Model.series` over the 600-epoch fixture **0.886 ms/epoch, so a day at 1 Hz is 1 min 16 s, not 2.1 hours**. **The cost is mostly a property of the BATCH, not of the toolbox**: fitting cost(K) = a + bK to the standard rows gives a ~= 15.2 ms of setup per CALL against b ~= 0.65 ms of actual per-epoch work (fast: 4.0 and 0.037). The 88 ms was a K=1 measurement and therefore almost entirely setup. The presets separate by 4.0x at K=1, 10.5x at K=100, and 17.6x on the marginal cost b -- which is the number that describes the physics, and where `geom.selfShadow` really does dominate. For the one caller that cannot batch, an integrator's inner loop, `app.DynamicsSource.accel` evaluates the whole standard chain for one state in **0.544 ms**, 29x cheaper than `Model.accel` at K=1. **A first attempt at this re-measurement was itself wrong** and is recorded in D.9: one warm-up call is not warm, and an argument built on the artefact ("the presets differ by 4 % at K=1") had to be withdrawn | closed | LP-R | — |
| O-27 | **Closed 2026-09-15 (LP-C WP-C11):** `env.wind("hwm14")` evaluates NRL's HWM14 + DWM07. Original entry: `env.wind` model `hwm14`: the Horizontal Wind Model 14 (Drob et al. 2015) source is not in hand, so the wind branch refuses. Quiet-time RMS is 40-80 m/s, which is 1 % of the along-track relative velocity and the dominant term in the cross-track drag, so this bounds what any cross-track drag claim can mean | closed | owner supplies HWM14, or LP-3b | — |
| O-22 | **CLOSED 2026-09-17, and it had been closed for eight days without anyone noticing.** The owner chose MIT on 2026-09-09 and `LICENSE` is in the repository, naming Matthias Weigelt; its own text excluded `mex/nrlmsise00/`, and that remainder was carried on as O-34. **O-34 closed on 2026-09-09 by deleting that directory** -- the source is fetched from its maintainer now and vendored nowhere -- so the exclusion covers nothing that exists. A remainder tracked in another row stops being visible the moment that row closes; the Status column added on 2026-09-17 is what surfaced it | closed | — | — |
| O-34 | **CLOSED 2026-09-09: the source is fetched from its maintainer, not redistributed here.** `mex/nrlmsise00/` is removed and gitignored; `buildMex` resolves `'MsisSource'` -> `NGF_NRLMSISE00` -> `mex/nrlmsise00` and REFUSES rather than fetching; **`setupNonGravForces` does the clone**, at the pinned commit `bc9a2fe`, only when the source is missing, and prints a line when it does. That division is the point: an installer may reach the network, a build step may not. The first attempt put the clone in a CI shell step and the self-hosted Windows runner has no bash (`bash: command not found`, 2026-09-09) — the right answer was not a different shell but the observation that setup already runs first from `ciGate`, is what a human runs anyway, and can own the URL and the pin as a single authority (R1); a YAML step would have been a second place to change them. **Two findings shaped it.** (i) The vendored copy could NOT be reproduced from upstream: its `nrlmsise-00.c` digest matches no commit in the maintainer's entire twelve-commit history, so it came from a local tree and nobody but its holder could obtain those bytes — which is what made pinning to it untenable for a public repository. (ii) The published oracle is unchanged: section 7 of the `DOCUMENTATION`, the 17-case table `tNgfEnv` parses, is byte-identical between the old copy and upstream, so the move does not disturb a single reference value. The owner ruled on 2026-09-09 to track upstream and take its bug fixes (Argon densities, a memory corruption, gtd7 at 32.5 km). Original entry: the right to redistribute `mex/nrlmsise00/` is unsettled, and MIT does not settle it. Read 2026-09-09: the only legal text in that package is section 1 of its `DOCUMENTATION`, a warranty disclaimer plus a request to send patches to the maintainer (Dominik Brodowski). It grants no redistribution right and none is implied by its absence. NRLMSISE-00 originates with a US government laboratory, whose work is often public domain, but *often* is not a legal finding. `THIRD-PARTY-NOTICES.md` records four ways out and recommends the cheapest: **remove the vendored source and have `buildMex` point at a user-supplied copy**, which costs a download and depends on nobody answering an email. `env.density` already refuses rather than degrades without the MEX, so only the fidelity tier is affected | closed | owner | — |
| O-23 | **CLOSED 2026-09-17, and it had been done since 2026-09-15 without the row being told.** `geom.MacroModel`'s unread per-panel `alpha` was removed in LP-C WP-C2 -- the file's own footer records it and the only mention left is a note saying it was required until that date. The row still described the toolbox as it had been, which is what a records sweep is for: a status column makes a table CONSISTENT, not TRUE. | closed | — | — |
| O-30 | **Closed 2026-09-08 the same day it was opened:** the MCP bridge never answered, but the CI `matlab` job on the self-hosted runner did. Predicted 270, ran 270, 269 passed; the one failure was a guessed absolute speed bound, decomposed in the same session. 0.6.1 is ☑ | closed | — | — |
| O-31 | **CLOSED 2026-09-18 with O-29: the question it asks no longer has a subject.** It asked for the kernel's share of `Model.accel`'s 88 ms, and there is no 88 ms -- WP-C3 took the standard preset to 0.798 ms per epoch at K=100, and most of what the old figure recorded was per-call setup rather than kernel work (see O-29). The ruling of 2026-09-08 (P-11) had already made this informational: both paths are kept by design, so no ratio between them decides anything. The numbers exist anyway, measured 2026-09-18: `kernel.singleEpoch` compiled 0.456 ms and the native whole-chain evaluator 0.544 ms, against the MATLAB standard preset's 0.798 ms at K=100 and 15.9 ms at K=1 | closed | LP-R | — |
| O-33 | **VOID 2026-09-09 with the owner's ruling that lsqvce is overkill.** WP-19 carries the toolbox's own solver, nothing external is called, and there is nothing for CI to check out. Original entry: WP-19 is blocked by the runner, not by design. `lsqvce` is checked out by neither CI job: the `matlab` job clones `referenceSystem` and `gracefoReading` only, with an organisation read token scoped to this repository. A `calibrate` test would therefore ERROR on CI rather than skip, which reddens the gate for an environment reason. Two ways out and the choice is the owner's: extend the runner checkout and the token scope to `astrobayer/lsqvce`, or accept that WP-19 is verified on the host and exempt from the CI suite — the second weakens the gate and should be a last resort. Nothing else in LP-7 depends on it. (Noticed 2026-09-09 while reading `dependencyReport`, whose lsqvce minimum still said 1.0 — a version whose surface no longer exists; corrected to 5.3 in the same commit) | void | owner | — |
| O-32 | **ANSWERED 2026-09-08 by measurement, sign ratio -1.0**: with `R = rotationFromQuaternion(q)` (body→GCRS), `omega = -referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference")`. The minus is not optional; it is written into the help of `sensor.Electrostatic` so nobody derives it twice. Original entry: the sign convention linking the two units is not settled by reading. That unit takes R the other way round (`Rdot = -Omega R`, i.e. R is navigation→body) and carries a configurable `.sense`, while `rotationFromQuaternion` on SCA1B returns body→GCRS, and the sensor's contract is fixed by `skew(omega) = R' Rdot`. The centrifugal term is even in omega and hides the difference; only the Euler term would be wrong, by exactly -1. `tests/tNgfSensor/ElectrostaticReference` therefore asserts the MAGNITUDES agree and RECORDS the measured sign ratio in its diagnostic instead of asserting a guess. Done | closed | — | — |
| O-37 | **A magnetorquer current is not yet a torque.** `readMAG1B` returns the six branch currents and their net per axis, and stops there: turning amperes into a magnetic moment needs the coil's turns-area product and its mounting, and neither is in the file nor in the handbook section this project has read. The torque itself is `m x B`, which is trivial once `m` exists. The conversion belongs to the magnetic-torque model (WP-23/24) and needs one number per axis from the GRACE-FO handbook or the manufacturer | open | owner or WP-23 | before a magnetic torque is modelled |
| O-36 | **`.mltbx` packaging is deprecated: the toolbox is published via git only** (owner ruling, 2026-09-09). `tools/packageToolbox.m` is removed rather than left unused, because a packaging script nobody runs is a second release path that will be believed. Nothing else in the repository depended on it; `setupNonGravForces` is the installer either way, and it is what the README tells a reader to run | closed | — | — |
| O-35 | **CLOSED 2026-10-18: answered by the owner 2026-09-09 and implemented the same day (`setupNonGravForces` 'FetchData'); its re-entry condition, the first `.mltbx`, cannot arise since packaging was deprecated (O-36).** Original entry: **ANSWERED by the owner 2026-09-09: data is optionally downloaded during setup, and otherwise the user supplies it.** `setupNonGravForces` gains `'FetchData'` with the same three settings as `'FetchMsis'`, fetching only when missing. **Corrected the same day: GRACE-FO Level-1B IS anonymously downloadable, and my probe was wrong.** I looked at `isdc.gfz-potsdam.de` and a guessed `isdc.gfz.de/grace-fo-isdc/level-1b/` path, got 404 and a login wall at PO.DAAC, and concluded no anonymous mirror existed. The owner supplied `https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/`; verified the same day, it lists one archive per day and the 2020-01-01 one holds every product under exactly this toolbox's fixture filenames. Setup now fetches it (147 MB, only when the members are missing, only the members kept). **The lesson is the one this project keeps relearning: two failed guesses at a URL are not a finding.** Original text: GRACE-FO Level-1B is not downloadable without an Earthdata account. Probed 2026-09-09 — IERS answered 200 with the fixed-format file, CelesTrak answered 503 once and 200 on the retry (which is why the fetch retries twice), and no anonymous GRACE-FO mirror was established. Setup names the source and the filenames it expects and downloads nothing. The L1B question is settled (see above). What stays open is whether `data/` should also leave the REPOSITORY as `mex/nrlmsise00` did. It was kept — the three providers publish redistribution terms where NRLMSISE-00 granted none, and removing it would make every test and demo depend on two networks staying up, one of which already answered 503 during the probe | closed | owner | — |
| O-24 | **VOID 2026-09-09**: the dependency was dropped before the design was written, so there is nothing left to design against. Original entry: discharged as a blocker 2026-09-08: the 5.3 surface was read (B.4) and `app.calibrate`'s shape follows from it — a kind-C batched declaration, priors declared in the columns `lsqvce.layout` announces, VCE by `Estimate="noise"`. What remains is the writing, which is WP-19 itself and not an obligation | void | LP-7 | — |
| O-19 | WP-34 GROOPS bridge: GROOPS is not installed on the owner's host (searched 2026-09-07) and its GPLv3 source may not be transcribed, so the black-box fixtures under `tests/fixtures/groops/` cannot be produced here. The owner runs `SimulateAccelerometer` for 2020-01-01 GRACE-C once and commits the CSV; the drag path is then a genuine cross-check, SRP/ERP/TRP an implementation comparison (both cite Wöske 2019)  **2026-10-26 (PR #107): TU Graz publishes GROOPS-processed `nonConservativeForces` (satellite frame, 10 s) for GRACE-A/B and C/D per day; excerpts for 2008-01-01 and 2020-01-01..03 are in `data/samples/itsg`.** **The series is MODELLED (owner, 2026-10-26)**: GROOPS-processed non-conservative forces for the real days, i.e. the black-box GROOPS model output this row asks for, published by TU Graz instead of produced here. **The comparison exists since 2026-10-26**: `tNgfModel/modelGroopsReference` reads TU Graz's nonConservativeForces, attitude and reduced-dynamic orbit of 2020-01-01 C with `gracefoReading.readGroopsInstrument` (minimum raised to 1.1), feeds their attitude and orbit into `Model.accel` and records correlation, RMS ratio and residual per axis. **Measured (CI 40406425492 and 40411254932)**, 600 x 10 s: correlation 0.983 / 1.000 / 0.968 (x, y, z), RMS ratio Model/GROOPS 1.50 / 1.18 / 1.20, magnitude ratio 1.119, RMS difference 0.39 of GROOPS's; radiation agrees at 0.945 / 0.999 / 0.994, the along-track difference is DRAG (model 2.07e-8 RMS against a GROOPS total of 1.75e-8, mean 1.67e-8 against 0.71e-8). The first run had TU Graz's quaternion unconjugated (it is bit-identical to SCA1B's; CLAUDE.md section 7). GROOPS itself follows ACT1B at 0.990 / 0.982 / 0.910 over the overlapping 30 min, asserted in the same test. **Beside it, the measurement**: `modelAct1bReference` gives 5.5 / 3.4 / 2.1 nm/s^2 RMS after one bias per day (Huckfeldt et al. 2024 eqs. 7-10), GROOPS against the same ACT1B 2.7 / 3.1 / 1.8 (stdlib computation), so this toolbox's along-track is 2x GROOPS's distance from the accelerometer. **Diagnosed 2026-10-27 (PR #109, `tNgfModel/modelDragDiagnosis`, CI 40447628923): the difference is a DRAG SCALE.** GROOPS = k_d drag + k_r radiation + k_t thermal of the standard preset gives k_d = 0.482 along-track, constant over the orbit (0.44-0.49 in four 25 min windows), and the fit takes the along-track RMS from 10.6 to 0.53 nm/s^2: GROOPS's drag is HALF of this toolbox's. No single option reaches it (msis21 0.915, Sentman alpha 1 0.905, constantCd 2.4 0.751 of the standard drag), and every option that lowers the drag also lowers the residual against ACT1B (5.57 -> 4.56 nm/s^2 with constantCd 2.4; msis21 + constantCd 2.4 together 0.685 of the drag and 3.98 nm/s^2). The measurement agrees in direction: with radiation and thermal held, ACT1B over the 3 x 30 min wants k_drag = 0.21 (RMS 5.57 -> 1.64 nm/s^2), weakly determined on 30 min arcs with per-day biases and ACT1B's own scale in it; step 2 (whole days) is what pins it. **Pinned 2026-10-27 (PR #110, `tools/dragScaleMonth.m`, run 40460797299 on the runner's DATAPOOL): all 31 days of 2020-01 GRACE-C at 10 s, THR1B firings cut (36-52 % of each day masked), one bias per day.** ACT1B wants a joint k_drag of 0.443 (per day 0.34-0.57, median 0.42) -- the same as GROOPS's 0.48 times ACT1B's own scale against GROOPS (0.92, stdlib on the excerpt): two independent references say this toolbox's along-track drag is about TWICE too large. Median daily bias-only RMS: x 6.36, y 2.84, z 3.60 nm/s^2 (Huckfeldt: 3-5 x, <= 2 y, z). Variants, drag ratio / median rms_x / joint k_drag: msis21 0.931 / 5.57 / 0.486; constantCd 2.4 0.742 / 4.41 / 0.555; msis21 + Cd 2.4 0.688 / 3.86 / 0.609 -- even the combination leaves a factor of 1.6. **Decomposed the same day (`tools/pyxcheck/gracefo_cda.py`, the independent kernel mirror): the areas are the handbook's (frontal plate 0.955 m^2), and Sentman at alpha 0.85 gives Cd*A = 3.4-3.6 m^2 unshadowed, i.e. Cd 3.5-3.8 on the frontal plate, of which ~30 % is the edge-on shear on the 16.7 m^2 of side, nadir and zenith plates; Huckfeldt's base case is Cd 2.4 on that plate (2.29 m^2).** The gas-surface model therefore accounts for 0.66-0.74 of the factor (0.742 measured through Model, with self- and flow-shadowing), and the rest, 0.44 / 0.74 = 0.6, sits in rho: the NRLMSISE-00 / NRLMSIS 2.1 density in the 2019-2020 solar minimum, the regime where the empirical models are known to run high. With Huckfeldt's own configuration (Cd 2.4, an MSIS density) the month gives 3.86 nm/s^2, inside their 3-5: the toolbox reproduces the published residual level, so no implementation defect is indicated. **ANSWERED 2026-10-27 (PR #111, DTM2020, month run 40470246020): the factor is the DENSITY, and the decomposition above was WRONG about the gas-surface model.** Same month, same method: `drag.density = "dtm2020"` (Sentman, alpha 0.85, the preset otherwise) gives a drag 0.529 of the preset's, a joint k_drag of **0.961** (per day 0.81-1.12) and a median daily along-track RMS of **1.79 nm/s^2** -- below Huckfeldt's 3-5. DTM2020 + constantCd 2.4 overshoots (k_drag 1.207, 1.93 nm/s^2). At 490 km on 2020-01-01 DTM2020 is 0.478 of NRLMSISE-00 and 0.522 of NRLMSIS 2.1 (`tNgfEnv/densityDtm2020Metamorphic`, CI 40470246057). So the MSIS density in the 2019-2020 minimum is ~2x high and Sentman's Cd 3.5-3.8 is what the data want once the density is right; Huckfeldt's Cd 2.4 was absorbing the density excess. The 'gas-surface model 0.66-0.74 of the factor' above was an inference from the ratio of two Cd values, not a measurement, and the measurement contradicts it. Open: whether DTM2020 should become a preset's density (owner; its licence is academic only), GRACE-D and another month (solar activity rising), and TU Graz's configuration for GROOPS's 0.48. **Owner, 2026-10-27: DTM2020's licence (academic, non-commercial only) rules it out as the way forward; it stays opt-in as the reference. The licence-clean route is an OBSERVED density scale: ESA's Swarm A/B/C POD-only densities (SW_OPER_DNSxPOD_2_, ESA open data policy), fetched through VirES HAPI by `setupNonGravForces('FetchPool')` (PR #112, first stage: data), a per-day rho_Swarm / rho_MSIS applied to the GRACE-FO drag (second stage), judged by the month workflow against DTM2020's k_drag 0.96.** The first fetch-data run (40494585676) brought 2020-01 for A, B and C (8640 records in the three fixture days each, A/C at 435-464 km, B at 502-533 km, some B records flagged anomalous). **MEASURED 2026-10-27 (month run 40500042274, 2020-01 C, 31 days): the observed density scale works without DTM2020.** rho_Swarm / rho_msis00 per day, median A 0.525, B 0.475, C 0.526, interpolated to GRACE-FO 0.478 (0.41..0.61 over the month); with the preset's drag scaled by it ACT1B wants a joint k_drag of **0.903** (Swarm B alone 0.916, A+C 0.834) and the median daily along-track RMS falls from 6.36 to **2.08 nm/s^2** (B 2.05, A+C 2.28) -- against DTM2020's 0.961 and 1.79, and Huckfeldt's 3-5. Mean heights A/C 445.8, B 513.1 km; GRACE-FO 502.7 km by |r| - 6371 km, which is the crude part (a few km against a geodetic height; referenceSystem's geodetic conversion is the refinement). Remaining ~10 % in k_drag: Swarm's local-time planes differ from GRACE-FO's, and the per-day factor is a daily mean. **Owner, 2026-10-27: it enters Model as an option** -- `drag.densityScale = "swarm"` with `env.densityScale` (PR #113), interpolating at GRACE-FO's geodetic height per epoch instead of the tool's |r| - 6371 km. **MEASURED 2026-10-27 with the option, both satellites (month runs 40514608768 C and 40514618274 D, 2020-01):** variant: joint k_drag / median daily rms_x [nm/s^2], C ; D -- preset 0.443 / 6.36 ; 0.479 / 6.21; dtm2020 0.961 / 1.79 ; 1.024 / 1.80; **`drag.densityScale swarm` 0.882 / 2.13 ; 0.937 / 2.11**; Swarm B scale 0.916 / 2.05 ; 0.976 / 2.08; Swarm A+C scale 0.834 / 2.28 ; 0.883 / 2.21. D: 27 days, 19-21 without records, 22 skipped (THR1B time_ref "G,R" after the safe mode); per-day k_drag with the option 0.79-1.07. The option's factor over the kept epochs is median 0.475 (0.41..0.61), GRACE-FO's mean geodetic height 506 km against the tool's 502.7 by |r| - 6371 km; its k_drag sits 0.02 below the tool's earlier height-interpolated 0.903 on C; the difference is the per-epoch height (the drag weights the low epochs), not decomposed further. The two satellites agree to 6 % in k_drag in every variant, so the scale is not a C artefact; D sits ~6 % higher throughout (its own drag ratio to ACT1B, not the density). **2026-10-27 (PR #114): a month of higher activity** -- `setupNonGravForces('PoolMonths', "yyyy-mm")` fetches a whole month into the pool (GRACE-C/D ACT1B, GNV1B, SCA1B, THR1B from ISDC's day archives; Swarm A/B/C; the complete SW-All and finals2000A, refetched when they end before the month), `fetch-data.yml` takes the months as a dispatch input, and `tools/dragScaleMonth.m` reads the pool's complete space weather and EOP and refuses a month they do not cover (it used the 2019-09..2020-03 excerpts, which serve 2020-01 only). Month chosen: 2024-01 (solar cycle 25 near maximum). First fetch (run 40569451735): all 31 L1B days, Swarm A/B/C, SW-All and finals2000A fetched -- ISDC answers again (O-51's GRACE-FO side). **MEASURED, 2024-01 GRACE-C (run 40574285164, 31 days, F10.7 ~165):** joint k_drag / median daily rms_x [nm/s^2] -- preset 0.964 / 20.95; msis21 1.060 / 20.59; constantCd 2.4 1.273 / 29.12; dtm2020 1.337 / 28.70; dtm2020 + Cd 2.4 1.790 / 49.31; **drag.densityScale swarm 0.978 / 19.50**; Swarm B 1.015 / 18.64; A+C 0.982 / 18.86. rho_Swarm/rho_msis00 median 0.955 (B) / 0.993 (A+C) / 0.986 (option), per day 0.77..1.15; heights A/C 488.6, B 521.7, GRACE-FO 494.5 km. **At high activity NRLMSISE-00 is right and DTM2020 is 0.77 of it, too LOW**: the 2020-01 factor of 0.48 is a solar-minimum error of MSIS, not a constant one, which is exactly what a driver model must carry. The Swarm option stays at k ~0.98 in both regimes (2020-01 0.882, 2024-01 0.978), DTM2020 moves from 0.961 to 1.337. The residual is ~10x larger than in 2020 (20 vs 2 nm/s^2), in proportion to the drag: within-day density variation that a daily factor cannot carry. THR1B masks only 0.3-22 % of a day in 2024 (2020: 36-52 %). **GRACE-D 2024-01 cannot be compared: every ACT1B_2024-01-DD_D_04.txt in ISDC's archives has no records** (run 40581217526; GNV1B, SCA1B and THR1B are present); the D side needs another accelerometer product for that period. **Owner, 2026-10-27: build a density correction that extrapolates in time in both directions**, in the recommended order: a driver regression of log f first (step 1), then a correction in exospheric temperature through the model itself, fitted to accelerometer densities of several missions and cycles with a secular trend, validated across cycles. Step 1 is `tools/densityDriverFit.m` (`.github/workflows/density-fit.yml`, on demand): per Swarm satellite, UTC day and orbit node the observed log f, regressed with `timeSeriesPreparation.regressDrivers` on Fbar (81-day mean ending the previous day), dF, ap, annual/semiannual and local-time harmonics, height and a trend, F10.7 against F30, nested; judged by yearly blocked cross-validation (`timeSeriesPreparation.crossval`, 27-day buffer) and by forward (fit before 2020, predict after) and backward prediction against the constant factor. Swarm 2014-01..2026-08 is fetched (2026-09 is not yet published: VirES answers 400 for a stop past the dataset, run 40574264167) by `setupNonGravForces('SwarmMonths')` (fetch-data input `swarm`). **STEP 1 MEASURED (density-fit runs 40597656794 and 40603989516, Swarm 2014-01..2026-08, 27 686 passes, 13 844 satellite-days):** the yearly median of rho_Swarm/rho_msis00 follows the solar cycle -- 0.84 (2014), 0.94 (2015), 0.82, 0.68, 0.59, **0.50 (2019)**, 0.52, 0.71, 0.80, **0.87 (2023)**, 0.77, 0.84, 0.86 (2026); A, B and C agree to ~1 %. 18 passes (0.07 %) are wild and take std(log f) from 0.25 to 0.65; trimmed at log 2 from the 27-day running median they are gone. On the trimmed days (RMS of held-out log f, i.e. relative error of a predicted daily factor; constant model in brackets): yearly blocked CV 0.217 (0.252) for F30 + season, backward (fit 2020-2026, predict 2014-2019) 0.220 (0.280), **forward (fit 2014-2019, predict 2020-2026) 0.287 (0.232) -- worse than the constant**. F30 beats F10.7 everywhere (CV 0.218 vs 0.224); the LT, height and seasonal terms add <1 %; a linear TREND ruins extrapolation (forward 0.81, backward 0.45) and is out. Full fit (F30, trimmed days, r2 0.42): Fbar +0.0058 +- 0.0011 per sfu, dF +0.0026 +- 0.0010, ap +0.0025 +- 0.0007, height -0.059 +- 0.019 per 100 km. Reading: ~20 % per-day scatter remains (Swarm noise plus MSIS's day-to-day error); the drivers carry the cycle; and the cycle-24 relation does not carry linearly into cycle 25 (f is not monotonic in Fbar at high activity: 0.94 in 2015 against 0.77 in 2024 at higher flux) -- the reason for steps 2 (a correction through the model's own temperature) and 3 (more missions and cycles). **NONLINEAR TEST (owner's choice, PR #115, run 40624842064, trimmed days; CV / forward / backward RMS of log f):** constant 0.252 / 0.232 / 0.280; linear F30 + season 0.217 / 0.287 / 0.220; **F10.7 + season + ((Fbar-120)/100)^2 0.192 / 0.291 / 0.192**; F10.7 + square + phase 0.197 / 0.252 / 0.195; F30 + square 0.179 / refused / 0.184 and F30 + square + phase 0.178 / refused / 0.189 (best CV; on 2014-2019 alone Fbar and its square are collinear, VIF > 10, so the forward fit is refused); log Fbar 0.203 / 0.254 / 0.206 (F30). The square is the real gain -- CV and backward 15 % below the linear model, 24-29 % below the constant -- and its coefficients (F10.7: Fbar +0.0039 +- 0.0004 per sfu, square -0.71 +- 0.11, both VIF 1.1) put the maximum of log f at Fbar ~148 sfu: MSIS's error is largest in the minimum, smallest near 150 sfu and grows again above, which the linear model could not express. The phase term adds nothing alone and little with the square. **No model beats the constant FORWARD**: a fit on 2014-2019 (Fbar up to ~150) extrapolates to 2024 (Fbar ~200) outside its range. Driver space, not time, is what must be covered -- the case for step 3 (cycle 23 densities from CHAMP/GRACE reach Fbar > 200) and for step 2 (a physical shape for the extrapolation). Open: which force models TU Graz's processing used, since a difference is only diagnostic once those are known| open **Step 3 begun (PR #116, 2026-10-28):** TU Delft's accelerometer densities of CHAMP, GRACE-A/B and GRACE-FO (ESA TOLEOS v02, CC BY 4.0) through the same VirES HAPI server: `setupNonGravForces('DensityMonths')` and FetchPool's reference-day months fetch them into `<DataPool>/TOLEOS/<dataset>/<yyyy>`; dataset names from viresclient's catalogue (`CH_OPER_DNS_ACC_2_`, `GR_OPER_DNS1ACC_2_`, `GR_OPER_DNS2ACC_2_`, `GF_OPER_DNS1ACC_2_`), verified by fetch-data run 40632322798, whose printed headers show Swarm DNSxPOD's nine columns and units at 10 s (CHAMP 2000-07-29..2010-09-04, GRACE-A 2002-04-04..2017-10-31, GRACE-FO 1 2018-05-29..2026-04-30; CHAMP's orbit mean is the fill 9.99e+32); GRACE-B left out as redundant; `tools/cutFixtures.py` cuts six hours of each into `data/samples/toleos`. `io.readSwarmDensity` reads the TOLEOS datasets (satellite CHAMP / GRACE-A|B / GRACE-C|D from the dataset name, values >= 1e30 NaN and not valid, an unknown dataset `dataSatellite`), four tests on the excerpts; `tools/densityDriverFit.m` takes 'Missions' (☑ CI 40636131057: 479/479, predicted 479). **The months (PR #117):** fetch-data run 40637542696 fetched 2000-07..2026-04 in 29 min except the six partial end months, refused with 400 because the request reached outside the dataset's startDate..stopDate; setup now clips the end months to the headers' spans. **The multi-mission fit (run 40651219543, 2000-08..2026-08):** 50376 passes over six satellites; the yearly median f is 0.57-0.59 in the 2008-09 minimum, 0.51-0.53 in 2019-20, 0.80-0.95 in the maxima. Trimmed days, CV / forward / backward RMS of log f: constant 0.234 / 0.228 / 0.248, F10.7 + season + square 0.178 / 0.181 / 0.186 -- the forward prediction across 2020 now beats the constant by 21 % (Swarm alone: 0.291 against 0.232), so covering driver space was the missing piece. log f peaks at Fbar ~150 sfu. F30 + square (best CV 0.171) is refused forward: its square is centred at 120 like F10.7's, collinear with Fbar on 2000-19; to be centred per proxy. **Centred per proxy (PR #118, run 40666688202):** centres F10.7 112.0 sfu, F30 75.7 sfu; F30's VIF 9.5 -> 1.4; F30 + season + square CV 0.172 / forward 0.181 / backward 0.180 against the constant's 0.234 / 0.228 / 0.248. Residual medians per satellite (log f): Swarm -0.003..-0.018, CHAMP +0.036, GRACE-A +0.045, GRACE-FO 1 -0.041 -- offsets of 3-5 %, a quarter of the per-day scatter; GRACE-FO 1 reads 4 % below Swarm in the same years. **The correction unit (PR #119, 2026-10-28, owner's request):** `env.densityDriverScale(sw, cls)` predicts f for any UTC day from F30 (81-day mean ending the day before, the day-before deviation), the day's Ap and the season, with the coefficient record of run 40676249479 (F30 + season + square, trimmed days, every digit; a0 -0.70150825262940808 absolute, centre 75.712339474163883 sfu, r2 0.517; maximum of log f at Fbar ~99 sfu of F30); 'Offset' adds a satellite's level in log f (GRACE-FO 1 -0.041). `drag.densityScale = "driver"` passes it to `Model` like the Swarm scale (S.kind must match the option). Oracle: the tool's own drivers and log f at MJD 58849 (-0.54833027703885062), computed from the full pool by a trailing movmean, against the unit's calendar window on the committed excerpts. Written without MATLAB, ◐, predicted gate 486 (479 + 7); first CI run 40682547930 485/486, the vectorization test red on `beta' * z` (a BLAS product sums in a batch-size-dependent order), now a column-wise sum; ☑ CI 40685698661: 486/486, predicted 486. **Judged by the month (PR #120, `tools/dragScaleMonth.m`, GRACE-C, joint k_drag / median daily rms_x nm/s^2):** 2020-01 (run 40708768247) preset 0.443 / 6.36, driver Offset 0 0.740 / 2.81 (drag ratio 0.601), Offset -0.041 0.771 / 2.63 (0.577), Swarm option 0.882 / 2.13 (0.493), dtm2020 0.961 / 1.79; 2024-01 (run 40708765369) preset 0.964 / 20.95, driver 1.072 / 20.34 (0.903) and 1.117 / 20.39 (0.867), Swarm option 0.978 / 19.50 (0.981), dtm2020 1.337 / 28.70. The driver scale, which needs no density observation of the month, removes about two thirds of the solar-minimum error in log k_drag and overshoots by 7-12 % at high activity; in both months its factor differs from the measured Swarm factor by ~0.2 in log f, the fit's own per-day scatter. The -0.041 offset helps in 2020 and hurts in 2024. **Step 2, the flux correction (PR #120, 2026-10-28):** `tools/densityDriverFit.m` finds per day the shift delta of NRLMSISE-00's F10.7 and F10.7A matching the observed sum and regresses it on the same drivers (fit run 40727237248): 13 % of trimmed days, all at solar minimum, lie below -60 sfu (MSIS cannot reach them); on the rest CV 0.178 against the constant's 0.194 (factor: 0.172 against 0.234); mission levels not closed (CHAMP -0.050, GRACE-A +0.034, GRACE-FO 1 -0.011 in log f). On ACT1B, GRACE-C, joint k_drag / median rms_x: 2020-01 flux 0.805 / 2.47, factor driver 0.740 / 2.81, Swarm 0.882 / 2.13; 2024-01 flux 1.095 / 20.19, driver 1.072 / 20.34, Swarm 0.978 / 19.50 -- the flux shift is as good as the factor or slightly better at GRACE-FO, worse in the fit statistics. **Owner's fallback built: `drag.densityScale = "combined"`** (`env.densityScale(..., 'Fallback', D)`; Swarm's factor on its days, D's elsewhere, per UTC day; ☑ CI 40747620816: 490/490, predicted 490). **GRACE-D 2020-01 with it (PR #121, run 40849007121, 27 days):** combined 0.937 / 2.11 nm/s^2 (joint k_drag / median rms_x), 0 of 120670 epochs from the driver, i.e. identical to the Swarm option; driver 0.802 / 2.59 (Offset -0.041: 0.836 / 2.45), flux driver 0.875 / 2.30, dtm2020 1.024 / 1.80, preset 0.479 / 6.21 -- the same ordering as on C, D ~6-8 % higher in k_drag. **Owner, 2026-10-29: the standard preset takes `drag.densityScale = "combined"`, strictly (no Model without 'DensityScale' evaluates it) and integrable: PR #122 gives the native path the daily table first.** | owner | before LP-R |
| O-25 | **CLOSED 2026-09-17 by comparing, which is the half of the choice that keeps the column honest.** `dependencyReport` measures `.belowMinimum` and still only reports (R5); `ciGate` refuses on it, as it does for `buildMex` and for the fetch report -- the gate is where a report becomes a verdict. The comparison is over the components the MINIMUM states, so "2.1" asks nothing about a patch and gracefoReading's "1.0.1" does. **Not every version is a number**: MICE answers `CSPICE_N0067`, which no dotted order can rank, so `.comparable` is reported and an unrankable version is never called too old -- inventing an order would have been a worse lie than the decorative column this replaces | closed | — | — |
| O-26 | **Largely discharged 2026-09-16, and what remains is recorded rather than fixed.** Tags now exist for `v0.1.0`, `v0.5.3`, `v0.6.0`, `v0.9.0` and `v0.10.0`. **Two things a later session must not be surprised by.** (i) `v.0.10.0` - with a stray dot - also exists, at the same commit `164ecdf`. It was a typo, and under the `~ALL` tag ruleset (deletion and non-fast-forward blocked, zero bypass actors) it cannot be removed: the one-way door this row warned about, walked through. It is harmless, points where `v0.10.0` points, and should simply be ignored. (ii) Only `v0.5.3` is ANNOTATED; the rest are lightweight, where `astrobayer/ci-actions/VERSIONING.md` asks for annotated tags plus a release. Those are immutable too, so the rule can only be honoured from the next tag onward. **And a sequencing lesson**: `v0.10.0` was cut at `164ecdf`, after which LP-D and the DE440 work merged WITHOUT a version bump, so `main` briefly served a changed `pilinski2010` constant under the same version string the tag names. Bumping to 0.11.0 closes that; the rule is that a tag marks a version and the next commit that changes behaviour opens the next one | open | owner | annotate from the next tag; cut it AFTER the bump, not before |
| O-40 | **`referenceSystem.downloadKernels` reports a digest that does not describe the file it wrote.** Its help says `.Digest` is the "SHA-256, once on disk"; for `de440s.bsp` it returned `addede59...164703`, while MATLAB's own `java.security.MessageDigest` over the written file, Python's `hashlib`, and the CI runner's independent computation all return `c1c7feea...260a49f2` for 32 726 016 bytes. Its `Mode "verify"` over the same CacheDir reports `Bytes NaN, Digest ""`, so it does not re-find the file it placed there either. This toolbox no longer quotes that number: `setupNonGravForces` measures the digest itself (LP-D WP-D4). **The lesson is more general than the bug**: a pin taken from the tool whose work it is meant to check is not a check. Not yet filed upstream - the sibling clone is on another session's branch and this project has not read that function's internals | open | owner / referenceSystem | file an upstream issue with the three measurements |
| O-41 | **CLOSED 2026-09-17 the same day it opened: the data is public and is now FETCHED, not vendored.** The owner named the source -- GROOPS publishes its whole data pool at `https://ftp.tugraz.at/pub/ITSG/groops/data/` -- so the redistribution question never had to be answered: nothing is committed. `setupNonGravForces` gains `'FetchGroops'` (default **never**, because two of the files are 1.25 GB EACH) and `'DataPool'`, and writes into `<DATAPOOL>/EarthRadiationPressure`. The pool is NAMED, never guessed: without the DATAPOOL environment variable the fetch reports and does nothing, because choosing a directory for 2.5 GB on a user's behalf is not a service (R5). No SHA-256 pin, deliberately: the publisher states it does not commit to keeping the files updated, so a digest measured once here would turn their next revision into a failure rather than a fact -- what is checked is the GROOPS magic and type string, read from the first bytes | closed | — | — |
| O-42 | **CLOSED 2026-09-17. Owner ruled WARN.** The preset half closed earlier the same day (`precise` takes March et al. 2021's fitted alpha = 0.85, WP-E7) and the citation was fixed in both kernel units. This closes the remainder: `kernel.accommodation` now reports `info.insideFit` for every model and raises `NonGravForces:accommodation:outsideFit` once per call when a computed model is evaluated outside the range it was fitted in -- alpha 0.85-1.00 for `pilinski20xx`, P_O >= 1e-5 Pa for `walker2014`; `fixed` and `sesam` never warn, the first because a stated value extrapolates nothing and the second because its source states no domain. Warn was chosen over refuse because extrapolating a fitted model knowingly is legitimate (R5). It WILL fire on GRACE-FO, which is the point | closed | LP-E | — |
| O-43 | **CLOSED 2026-09-18. (a) fixed; (b) WAS NOT A DEFECT, and the record claiming it was is corrected here.** (a) `ciGate` snapshots `pathdef.m` before installing dependencies and restores it afterwards, so the SESSION path keeps what the suite needs and the MACHINE keeps what it had; verified by hashing `pathdef.m` across a full gate run, byte-identical, with the gate printing that it had restored it -- the half that matters, because an unchanged hash alone is also what a restore that never fired looks like. The re-entry condition this row carried ("pass 'Save', false where the dependency supports it") could never have worked: `setupGraceFOReading()` takes no arguments and calls `savepath` unconditionally. (b) **This row claimed `setupReferenceSystem` warns `mexBuildFailed` "whether the binary is stale or merely locked-and-current", and that a consumer cannot tell them apart. MEASURED 2026-09-18: false on both counts.** That setup decides by DIGEST, not by timestamp -- it compares `mexGateway.sources` against `mexSourceDigest` and rebuilds only on a difference -- and on this host the digests genuinely differed (recorded 109793179 against 1898687363 in the isolated clone and 1010033319 in the sibling), so the binary really WAS stale and the warning was telling the truth. And the rebuild it attempts SUCCEEDS in a cold process: 20.6 s, 257 sources, manifest updated, digest matching afterwards. It failed during this session's gate runs because the gateway was already mapped, which this session's own concurrent MATLAB processes explain. **No upstream issue was filed, because there is nothing upstream to fix.** The rs_main gateway has been rebuilt and is no longer stale | closed | LP-R | — |
| O-44 | **CLOSED 2026-09-18 on the ruling, with the justification corrected: it was false.** This row said `ceresClimatology` has no native evaluator because a 64800-cell grid "is not table material for an integrator's inner loop". **The native path already carries THREE EBAF grids of exactly that size**, plus their axes and month nodes, for `erp.backend = "ceres"` -- 3.1 MB against the 2.1 MB this backend would need, and `ngfTiles` already does the bilinear lookup. Memory was never the obstacle and the port is about forty lines. The REAL reason to refuse is that it would buy no physics: `io.importCERES` of EBAF 4.2.1 is also 180x360 at one degree and also a measured exitance (verified 2026-09-18), so over an integration span the two differ in their source and their reference level (20 km against 30 km), not in what they can express. What `ceresClimatology` uniquely offers is an epoch range of 1900-2101 and provenance shared with GROOPS, and both are post-processing concerns. **An integrator uses `ceres` where the months are in hand and `knocke` where they are not.** The refusal message and the user guide now say that instead of the memory claim. Reopen in one sentence if an integrated arc ever needs an epoch outside the EBAF months held -- that is the one capability the port would add | closed | LP-R | — |
| O-39 | **CLOSED 2026-09-16, the same day it opened: the owner widened the secret's repository access and the runner checks odeSuite out.** The wall had been hiding two real defects, both invisible on this host and both fixed in the same round: (i) `buildMex`'s walk to odeSuite's `csrc` overshot by one directory, and since CI sets no `NGF_GJTOOLBOX` the `rhs` target would have reported unbuildable and `ciGate` would have gone red blaming an absent gjToolbox that was sitting right beside it; (ii) `ciGate` never built gjToolbox's own native stepping kernel, which is not committed, which a fresh checkout therefore lacks, and which `tNgfIntegration` propagates through. Each was masked by something this host had and the runner does not — an environment variable and a hand-built MEX. **The access wall was worth hitting: it is the only reason either defect was found before a release.** Original entry: **an ACCESS fact, measured twice, not a guess.** Since LP-C the `matlab` job checks the sibling out, because `app.DynamicsSource` subclasses `gjToolbox.DynamicsSource` and cannot load without it. `REFERENCE_READ_TOKEN` syncs `referenceSystem` and `gracefoReading` in that same job, seconds earlier, and then: **404** at the REST call that resolves the default branch (Metadata:read), and — with `ref: main` given so that call is skipped — **403** on the git fetch itself (Contents:read). Both permissions are refused, so the secret's repository scope does not include odeSuite. **The fix is one setting and it is the owner's**: add `astrobayer/odeSuite` to that organisation secret's repository access. The alternative, `continue-on-error` on the checkout, is refused here: the job would then run without gjToolbox and report green about a gate that never loaded the unit it exists to cover. Until then LP-C is ☑ on Tier A and unseen by the runner | closed | owner | — |
| O-38 | **CLOSED 2026-09-16 by removing the second authority, not by rebuilding it.** The investigation found the obligation had understated the problem: the committed PDF was the one from the initial commit, its title page said **version 0.1.0** three weeks and nine minor versions later, and the `.tex` it came from carried **eight statements the shipped toolbox contradicts** — `sesam`/`langmuir` accommodation (never implemented, O-16), DTM2020 and JB2008 as density backends (removed in WP-C2), `erp.anisotropy=adm` (removed, F4), an icosphere LUT (a cube map since WP-C3), complex-step verification (measured impossible and removed), and `lsqvce` as the estimator (dropped 2026-09-09). Rebuilding the PDF would have rendered every one of them faithfully. **The fix is the mechanism, not the file:** (i) `tools/buildDocs.py` generates `doc/src/generated.tex` — the version, the date and every enumerated option list, the last read from `Options.validateValue`'s own `sets` struct — and the document `\input`s it, so it cannot print a version or offer an option the toolbox would refuse; (ii) `auditPackage.py` refuses a stale `generated.tex` by CONTENT, refuses a committed PDF under `doc/`, and refuses a version literal in the `.tex`; (iii) a new CI job `theory` compiles the document on every run with `-halt-on-error`, fails on an undefined reference or control sequence (LaTeX only warns and still writes a PDF), and uploads the result as the `theoryDocument` artefact. The eight stale statements are corrected and the document gained the section it most lacked, `sec:flow`, on why occlusion is not facing and a gas is not a beam. Original entry: the theory PDF is one revision behind its own source; this host has no LaTeX and no instrument inspects the PDF | closed | — | — |
| O-11 | **CLOSED 2026-09-09 (WP-21).** `tools/buildDocs.py` renders all seven pages and `helptoc.xml`, each derivable fact from its single authority: the catalogue for the reference, the sensor table and the delivered/stub counts; `fail.m`'s help block for the mnemonic taxonomy; `toolboxVersion.m` for every footer. `auditPackage.py` refuses a stale page by CONTENT, and was fault-injected to prove it. The obligation had already come due: `Sensors.html` credited `Electrostatic` with a "CRN filter" — the one thing that class explicitly does NOT build — and listed all three sensors as future work months after two had shipped. A page nobody generates is a page that lies eventually | closed | — | — |
| O-49 | **CLOSED 2026-10-15 by release 0.15.0.** Version bumped (MINOR: N-03 added the time-varying scale and three units -- compareCalibrations, perRevolution, attributeRevolutions); `tools/benchmark.m` re-run on Tier A over fresh dependency clones and `tools/benchmark_reference.csv` committed at 0.15.0: standard K=1 17.0 ms, K=100 0.836 ms; fast K=1 5.15 ms, K=100 0.096 ms; compiled 0.521 ms, native 0.505 ms -- within a few per cent of 0.13.0's quiet run, confirming that 0.14.0's 5-30 % higher table was the loaded host. Gate 405/405, predicted 405. Original entry: **Release 0.15.0 (N-03, 2026-10-15).** N-03a-c are on `main` (after #85) at 0.14.0: the time-varying scale, model comparison, the per-revolution fit and its attribution -- added functionality, so MINOR. As O-48: bump `toolboxVersion`, every `Contents.m`, README, the HANDOVER title (audit check 5 compares it), re-run `tools/benchmark.m` on Tier A and commit `tools/benchmark_reference.csv` | closed | owner | — |
| O-48 | **CLOSED 2026-10-15 by release 0.13.0.** `tools/benchmark.m` was re-run on Tier A (R2026a Update 4, over a dependency root of fresh clones at origin/main, median of 5 after five warm-ups) and `tools/benchmark_reference.csv` committed at 0.13.0. Against the 0.12.0 table the MARGINAL per-epoch cost is unchanged (standard 0.646 -> 0.667 ms, fast 0.0367 -> 0.0349 ms, two-point fit at K=1 and K=100) while the fixed cost per CALL rose by about 2.5-2.8 ms in both presets (fast K=1: 4.00 -> 6.54 ms; standard 15.9 -> 18.7); the compiled paths did not move (0.456 -> 0.457, 0.544 -> 0.518). One run, so not yet attributed: the 0.12.0 table was taken against the referenceSystem of 2026-09 and this one against 2.3. The `Model.series` day figure (1 min 16 s) is produced by no instrument and was not regenerated. Original entry: **Release 0.13.0 (N-02, 2026-10-15).** N-02a-c are on `main` at 0.12.0; the version bump (toolboxVersion, every Contents.m, README) waits only for `tools/benchmark.m` to be re-run on Tier A and `tools/benchmark_reference.csv` committed, because audit check 13 refuses a reference measured against another version. Then O-47 (delete the qc forwarders) follows in the next PR | closed | owner | — |
| O-47 | **CLOSED 2026-10-15 by release 0.14.0.** The six forwarders, `+qc/Contents.m`, `tests/tNgfQc*.m` (39 tests), the catalogue rows, the hushes, `eventTable`'s exemptions, the L7 layer lines and the two orphaned oracles are deleted; callers in help text name `timeSeriesPreparation.qc`; 54 units. Original entry: **Delete the `qc` forwarders (N-01, 2026-10-14).** The six `NonGravForces.qc` units forward to timeSeriesPreparation for ONE release and then go, with `tests/tNgfQc*.m`, their catalogue rows, the six `deprecated` hushes in `NgfTestCase`, and `CLAUDE.md`'s L7 line (CONCEPT 10.1 N-01). `io.readTHR1B` already calls the target directly, so nothing in the toolbox depends on the forwarders | closed | session | — |
| O-46 | **Classic GRACE mass (N-07, 2026-09-22).** R-5 Sec. I.5.3.3 gives only the 487.2 kg LAUNCH mass and points to MAS1A/MAS1B for the mission value, and for classic GRACE that product carries the SPACECRAFT mass (RL02 `MAS1B_2008-01-01_A_02`: 476.85 kg, measured in gracefoReading's fixture), the opposite of GRACE-FO's tank mass. `io.massModel` serves GRACE-FO only and refuses anything at or below 500 kg, so `geom.grace` carries the launch mass and a classic user must pass `'Mass'` to `Model`. R-5 also prints no inertia tensor and no CoG -> proof-mass offset; both are NaN in `geom.grace` | open | N-08 | when N-08 builds the classic presets (#65) |
| O-45 | **Accelerometer calibration and the shared preparation toolbox (owner, 2026-09-18).** Design and decisions: the concept, `astrobayer/timeSeriesPreparation` `docs/CONCEPT.md` (Rev 0.12; moved there from gracefoReading at X-00a, timeSeriesPreparation#2), section 11 Q2, Q7, Q10, Q14, Q15, Q16; order and procedure in `docs/GUIDE.md` of the same repository; not restated here. **Consequence for this repository:** `app.calibrate` stays and is upgraded in place; `qc` leaves for timeSeriesPreparation; lsqvce re-enters the dependency table (B.4's own re-entry condition). Work, in the guide's order, one issue each: N-04 #61 → N-01 #58 → N-02 #59 (+ N-06 #63) → N-03 #60 → N-11 #68 → N-12 #69 → N-09 #66 → N-10 #67 → N-05 #62; classic GRACE N-07 #64 → N-08 #65. O-1 and O-18 (`qc` speed) move with `qc`: X-05 ports it unchanged and X-16 makes it fast there, so neither is worked on here. Until N-01 and N-02 land, do not extend `qc` or `app.calibrate`'s solver: that work would be deleted | open | owner | timeSeriesPreparation 0.1 (N-04), 0.2 (N-01), 1.0 (N-02), 1.1 (N-09) |
---

## Part I — Decision log

| Date | Decision | Reason |
|---|---|---|
| 2026-10-27 | **DTM2020 (operational: F10.7, Kp) enters `Options` as `drag.density = "dtm2020"`, opt-in, in no preset; supersedes the ruling of 2026-09-15 that took it out** (owner asked for it after the drag diagnosis) | O-19's month put ~0.6 of the along-track factor in the MSIS density; DTM2020 is fitted to TU Delft's CHAMP, GOCE and Swarm densities, which were derived with high-fidelity satellite models, so it is the independent density the question needs. CNES's licence (academic, non-commercial only; no modification, no redistribution) is met the way NRLMSIS 2.1's is: fetched by setup at a pinned commit and SHA-256, compiled unmodified, licence copied beside the binary, nothing committed. The operational version, not the research one: its drivers are in the CelesTrak file in hand, while the research version needs Hp60 (not in hand) and ships a reader that cannot be called unedited (O-52) |
| 2026-10-26 | **The model is compared with ACT1B the way Huckfeldt et al. (2024, ASR 73, Sec. 3.5) do: one bias per day and axis, scale 1, thruster firings cut; a common scale with per-day biases is recorded beside it, never instead** (owner asked for the comparison after the GROOPS one passed) | A bias-only fit keeps model error out of the calibration parameters, so its RMS is a statement about the MODEL and is comparable with the published 3-5 nm/s^2 along-track. THR1B is committed for 2020-01-01 only: that day is masked over [t - 0.5 s, t + onTime + 7 s] (their Fig. 3) widened by the CRN half span 70.3 s, the other two days are not, and day 1 unmasked is recorded to show what that costs (2 % in x, 27 % in y) |
| 2026-10-26 | **O-19's comparison reads TU Graz's GROOPS files through `gracefoReading.readGroopsInstrument`, not a reader of this toolbox; the gracefoReading minimum rises from 1.0.1 to 1.1** (owner asked for that reader) | A.2: this toolbox implements nothing it can call, and gracefoReading reads the real GROOPS dialect since 1.1.0. 1.1 contains the 1.0.1 fix of O-12, so the PATCH exception to VERSIONING.md ends rather than moves |
| 2026-10-26 | **All observation data lives on the CI runner, and setup collects it** (owner): `setupNonGravForces('FetchPool', 'PoolDays')` fetches Kyoto final Dst, CLS F30 and TU Graz's attitude, nonConservativeForces and reducedDynamicOrbit into DATAPOOL (`SpaceWeather/`, `ITSG/operational/`); `tools/cutFixtures.py` (no network) cuts the committed excerpts from the pool; `.github/workflows/fetch-data.yml` runs both on the self-hosted runner into `READING_DATA_POOL` and only the commit job runs on GitHub's machine. ciGate passes `'FetchPool', "never"`. Classic GRACE L1B is parked (O-51) | The development container cannot reach any of these hosts and GitHub-hosted runners cannot reach GFZ. A first version fetched with a Python script on a hosted runner, which made it a second network authority beside setup (R1); setup owns every URL now. First run 40381190966: all 30 products fetched, the 24 known ones with unchanged SHA-256s, every committed excerpt reproduced byte for byte |
| 2026-10-26 | **Kyoto FINAL Dst and CLS F30 are the sources for O-9; OMNI2 and Nobeyama's raw file are not** | Kyoto is the primary source of Dst and has final values for both fixture months (provisional/realtime answer 403 for them); OMNI only copies it. Nobeyama's own daily file has NaN at 1 GHz on 2020-01-01..03; CLS interpolates and flags those days, and DTM2020 is driven by the CLS series |
| 2026-10-26 | **TU Graz's `nonConservativeForces` is MODELLED** (owner, later the same day; supersedes the row's first ruling, kept below as it was) | The owner settled what TU Graz's README leaves open. It is therefore an implementation comparison for O-19 (GROOPS against this toolbox on the real day), not a measurement |
| 2026-10-26 | **(superseded the same day) TU Graz's `nonConservativeForces` is NOT an O-19 oracle yet** | Their README gives the frame (satellite reference frame, 10 s) but not whether the series is the calibrated accelerometer or a model; the values (~1e-8 m/s^2, no bias) rule out raw ACT1B and nothing else. Settle it with TU Graz before any comparison |
| 2026-10-20 | **ACH1B is an oracle, not a baseline: regressionBaseline runs it for GRACE-D in trial mode only, its arcs labelled ACT1B for timeSeriesPreparation's preset; the measured comparison is recorded, not frozen** (N-10b, ratified by the owner 2026-10-20) | the frozen baseline must stay one product to compare releases; the label is the smallest honest workaround until timeSeriesPreparation's Config knows ACH1B |
| 2026-10-20 | **A temperature-dependent scale S(T) = S + S1 (T - T0): S1 on the estimated entries only, one reference T0 over ALL arcs and S1 common in multi-day fits, a high VIF of its columns a warning not a refusal, refused beside a time-varying scale** (N-11b, ratified by the owner 2026-10-20) | a per-day reference would delete the day-to-day temperature differences that are S1's only leverage on a thermally controlled instrument; S1 is a physical coupling like the scale; time and temperature dependence are not separable on one arc |
| 2026-10-20 | **A measured series is never interpolated across a gap longer than 10x the 90th percentile of its sample spacing** (N-11b, also applies to N-11a's regressors, ratified by the owner 2026-10-20) | a 20 min gap interpolated moved S1 by 3.7 %; the median fails on burst-sampled MAG1B and AHK1B |
| 2026-10-19 | **The day bootstrap of a multi-day calibration holds the full fit's noise and robust weights fixed in every replicate (lsqvce.solveContributions); lsqvce minimum 6.35** (ratified by the owner 2026-10-20) | under the defaults (white, no robust) nothing changes -- identical sigmas to ten digits on real data -- and the bootstrap became ~4x cheaper per calibrate; re-estimating weights per replicate would need a full refit each time again |
| 2026-10-17 | **The standard preset reads the gas-stream shadow from geom.flowShadow's table (Options geom.flowShadow "table"); precise stays exact; the native path reads the same table** (owner, 2026-10-16: "use default as standard") | it moves the drag by at most 6.7e-13 m/s^2 on ten real days and halves Model.accel's cost; precise keeps the computation it tabulates |
| 2026-10-16 | **The regression baseline is the `standard` preset, run at a release on Tier A (~4 h), compared at k = 3 sigma (scale: bootstrap over days; bias: formal); the audit refuses a version it was not verified at** (N-10a, ratified by the owner 2026-10-16) | `fast`'s guessed density gives a meaningless scale (0.32); the month does not fit a 45-minute CI job; the audit is what makes a release-time instrument actually run |
| 2026-10-15 | **Measured regressors are a general 'Regressors' Arc, interpolated linearly, demeaned per arc, filtered like the bias; per-day columns with one arc, COMMON coefficients over several; apply refuses a calibration with regressor terms when the series are not passed** (N-11a, ratified by the owner 2026-10-16) | one mechanism serves temperatures, torquer currents and any later housekeeping; linear does not overshoot between housekeeping samples; dropping a fitted term silently would apply a different calibration |
| 2026-10-15 | **A multi-arc calibration's common scale takes its uncertainty from a bootstrap over DAYS (block length 1, each replicate refitted jointly), and warns below 10 days; the per-day biases carry lsqvce's formal local covariance** (N-09; RATIFIED by the owner 2026-10-15) | days are independent where the samples of one day are not, so whole days are the honest resampling unit; a per-day bias is not defined in a replicate that drew its day twice or not at all |
| 2026-10-15 | **The joint fit supports a constant or linear bias per day, a constant common scale and no jumps; everything else is refused, with 'ScaleMode' "perArc" as the route** (N-09; RATIFIED by the owner 2026-10-15) | a 'day' local block must have one size on every day; a common curve over days is undefined here; per-day jumps change the block size |
| 2026-10-15 | **`setupNonGravForces`'s DE440 fetch reads `referenceSystem.downloadKernels`' manifest before it hashes, and `localSha256` refuses a file it cannot open before it registers `fclose`.** A FAILED row, or a status with no file written, is recorded as `failed: download: <reason>`; a present file that cannot be opened is `failed: cannot open ...` instead of an uncaught error. `ciGate`'s `dataMissing` message now says "download or content check" | CI run 38733381139 (PR #85): NAIF closed the connection, downloadKernels warned and wrote nothing, and hashing the missing file raised "Invalid file identifier" from the onCleanup `fclose(-1)`, so the gate went red under the wrong cause. FAULT-INJECTED with a scratch shim of `downloadKernels` and the kernel moved aside: a FAILED manifest gives `de440: failed: download: injected: NAIF closed the connection ...` at setup AND as ciGate's `NonGravForces:ciGate:dataMissing` message; a "fetched" row with no file gives `failed: download: ... wrote no de440s.bsp`; a kernel held under an exclusive lock gives `failed: cannot open ... to hash it`; with the shim removed, the same call reports `present`. Gate over fresh clones: 391/391, predicted 391, codeIssues 0; after merging N-03c from main, 405/405, predicted 405 |
| 2026-10-15 | **The attribution's verdict counts a driver at 3 bootstrap sigma by default, not 2; a constant or missing driver is dropped and named; collinearity is the verdict "unresolved", not an error** (N-03c; RATIFIED by the owner 2026-10-15) | measured: the block bootstrap on 60 revolutions gave a null driver 0.54 of its OLS sigma, so at 2 sigma both v3 cases read "mixed"; at 3 sigma 12 of 12 were right. CONCEPT 6.2 makes the VIF refusal a verdict |
| 2026-10-15 | **Per-revolution windows are fixed in time from the first usable epoch, not tied to the ascending node; a window is fitted with a constant bias and a diagonal scale through fitCalibration, and a gated window is NaN with its reason** (N-03c; RATIFIED by the owner 2026-10-15) | v3 WP-16; one estimator applied per window; node-tied windows need the orbit in perRevolution, which the observation-only signature avoids |
| 2026-10-15 | **Model comparison ranks by blocked CV only; AIC/BIC are reported as IC_unreliable reference columns, never used to rank** (N-03b) | v3 WP-14; measured: AIC and BIC chose the overfit spline on 9 of 9 axes where CV chose the planted model on 9 of 9 |
| 2026-10-15 | **The CV buffer is at least the CRN window; the minimum-data gate scales by the training share on a fold; a local basis whose knots do not exceed the held-out window is refused** (N-03b; RATIFIED by the owner 2026-10-15) | a filtered row shares 140.7 s of raw samples with its neighbours; the gate judges data, not the CV design; a fold that cannot fit a candidate must not silently drop out of its score |
| 2026-10-15 | **A time-varying scale is identified against the bias by BASIS SIZE, not by a penalty: a knotted scale may not be finer than a knotted bias** (N-03a; RATIFIED by the owner 2026-10-15) | v3 WP-14 asked for an asymmetric smoothness prior with lambdas by blocked CV; a coarse basis states the same asymmetry without a weight to tune, and lsqvce priors stay the route if a penalty is ever needed |
| 2026-10-15 | **`isDriftSignificant` is a Wald test of "the curve is constant", and it WARNS on the formal white-noise covariance** (N-03a; RATIFIED by the owner 2026-10-15) | v3 named "drift amplitude vs band" without a criterion; the Wald test is that comparison made joint over the curve. Measured: the formal covariance flagged a planted-zero drift at p = 1e-18, the bootstrap one gave p = 0.72 |
| 2026-10-15 | **N-03 is three PRs under #60: a time-varying scale, b model comparison, c per-revolution and attribution** | the N-02 pattern; each part has its own oracle and the attribution needs the per-revolution fit first |
| 2026-10-15 | **O-47 ships as release 0.14.0 (MINOR), and L7 is left vacant rather than renumbered** | VERSIONING.md: `0.y.z` carries no stability promise, so a removal is the next MINOR; renumbering `app` to L7 would falsify every "L8 app" already written in records nobody compares |
| 2026-10-15 | **A.2's rule reads "NonGravForces computes forces and torques, and applies them; it implements nothing it can call"** (N-06; RATIFIED by the owner 2026-10-15) | the app layer calibrates against the forces (O-45); the old wording forbade that, the new one keeps what the old protected: no gravity, no integrator, no least squares inside |
| 2026-10-15 | **`app.calibrate` evaluates a Model at the ORBIT's epochs and interpolates onto the 10 Hz raw grid by pchip; jumps from run's QC are steps shared by the three axes, fitted only on request ('Jumps' true; measured: by default they were spurious)** (N-02c; RATIFIED by the owner 2026-10-15 -- the jump half directly, the orbit-epoch evaluation with pchip after the owner asked why and it was measured: the owner ruled "the measurement settles it". Its accuracy claim had been argued, not measured; MEASURED 2026-10-15 across a shadow entry against a 10 Hz evaluation, through calibrate's own route: 1.1e-10 m/s^2 raw, 9.3e-13 after the CRN filter, 3e-5 of the signal; a sample-and-hold fault reads 2.0e-2 and fails the 1e-3 bound (tNgfCalibrate/calibrateReferenceInterpolationPrecision)) | 10 Hz evaluation costs ten times as much for content above the CRN band; pchip does not overshoot an eclipse transition. A jump on one axis costs the others one parameter each, which is cheaper than a per-axis basis |
| 2026-10-15 | **Release 0.13.0 is its own PR, after N-02c** (O-48) | the audit refuses a version whose `tools/benchmark_reference.csv` was measured against another, and that file is written by MATLAB only |
| 2026-10-15 | **N-02b's estimator is its own unit, `app.fitCalibration(obs, ref)`; its noise default is WHITE with a moving-block bootstrap, not an estimated law** (taken in N-02b; RATIFIED by the owner 2026-10-15) | a separate unit lets N-02c wrap the numerics without re-opening them. The CRN-filtered residual is band-limited at 0.035 Hz, so a whitening law estimated from it (lsqvce's `Estimate "noise"` or an ARMA from `timeSeriesPreparation.noise`) is ill-conditioned above the band, while the pairs block bootstrap needs no law; both laws stay available as `'Noise'` options. The Calibration carries the bootstrap covariance |
| 2026-10-15 | **lsqvce re-enters the dependency table at minimum 6.22, optional for the toolbox** (N-02b) | CONCEPT Q2: the one estimator. 6.22 is the version read and gated when N-02b was written; MAJOR.MINOR per VERSIONING.md |
| 2026-10-14 | **N-02 is three PRs under #59 (a: truth + `app.Calibration`; b: estimator on `lsqvce.fit`; c: arc entry, run, attribution, release 0.13.0), and `app.calibrate`'s old signature is replaced outright in 0.13.0** (owner, chosen from options put) | each PR is reviewable and gated on its own; a deprecation shim would keep two calibration APIs alive for a release with nobody using the old one |
| 2026-10-14 | **`sensor.Electrostatic`'s default kernel is `timeSeriesPreparation.crnKernel("ACC")` on 10 Hz dynamics; a 1 Hz arc needs an explicit `[]`** (owner, N-04, chosen from three options put) | the concept row names the CRN kernel as the default; a rate-dependent default would degrade silently, an opt-in would not be a default. Supersedes the 2026-09-08 row "The CRN filter is an argument" for the default only |
| 2026-09-22 | **`geom.grace` carries only what R-5 prints; unknown mass properties are NaN, not zero or GRACE-FO's** (N-07) | the concept requires the macro model from the supplied document, never recalled (V-9); a NaN inertia makes the gravity-gradient torque refuse, a borrowed one would make it silently wrong |
| 2026-09-18 | **The concept (`astrobayer/timeSeriesPreparation` `docs/CONCEPT.md`) is the design authority for accelerometer calibration and preparation** (owner; decisions in its section 11, not restated here) | consequences for this repository are O-45 |
| 2026-09-16 | **Both Pilinski fits ship, as separate options** (owner): `pilinski2010` takes K = 7.50e-17 m^3/K from Pilinski, Argrow & Palo (2010) eq. 8, `pilinski2011` takes K = 4.98e-17 from Pilinski (2011) sect. 7.3 | the owner supplied the 2010 paper and it fits a DIFFERENT constant from the one shipped under its name. One option cannot carry two published values; naming each after the source that fits it lets a user reproduce either result, which is what a reference implementation is for |
| 2026-09-16 | **The papers close neither O-16 nor O-17, and that is recorded rather than stretched** | Pilinski 2010 is the Langmuir accommodation paper, not SESAM 2013: it contains no sticking law, no shape factor and no substrate model (the word SESAM does not occur). Walker, Mehta & Koller 2014 JSR 51(5) A32668 computes CLL by DSMC sampling of the scattering kernels and gives no closed form; the closed-form CLL model is the companion paper A32677 of the SAME issue. What they DO give is better than a partial implementation: published closed-form oracles and a quantified bound on the approximation this toolbox ships |
| 2026-09-15 | **LP-C opened; LP-R waits for it** (owner) | a review of the whole repository found fifteen defects (LP-C, Part G); the release is cut from a toolbox whose presets evaluate and whose integrator entry is usable |
| 2026-09-15 | **Self-shadowing is occlusion only; facing is the kernels'; gas occlusion is flux-weighted over the incident Maxwellian** (owner approved the design as proposed) | the table's "f = 0 on unlit plates" convention was a radiation-beam rule applied to a thermal gas, and it deleted the edge-on shear that Sentman's integrals carry on ~15 m² of GRACE-FO plates; the same convention zeroed the Earth heating of every panel facing away from the Sun (F1, F2) |
| 2026-09-15 | **`gjToolbox.DynamicsSource` is the integration contract, with a native C entry** (owner) | gjToolbox's interface exists and is frozen (read at odeSuite `origin/main`, `csrc/gj_native.h` ABI 1); the handle `Model.rhs` returned was never callable by it (F3) and is removed rather than kept beside the source |
| 2026-09-15 | **The precise preset is finished, not removed; data is fetched or its source named** (owner) | F5: a documented tier that cannot evaluate one epoch is a claim, not a feature |
| 2026-09-15 | **HWM14 is delivered as an opt-in; DTM2020, JB2008, CLL and Schamberg leave Options** (owner) | winds dominate cross-track drag and HWM14's source is open; DTM2020 needs an F30 source nobody holds, JB2008's code terms are unread, CLL has no closed flat-plate form. A value Options accepts must evaluate (R4) |
| 2026-09-15 | **CERES: EBAF-TOA monthly now (obs4MIPs, anonymous), SYN1deg-3Hour when the owner supplies a file** (owner) | EBAF is fetchable and testable today; SYN1deg resolves the diurnal cycle but is behind an Earthdata login the assistant does not hold; its variable names are read from a real file, never assumed (O-3) |
| 2026-09-15 | **Geodetic conversion and quaternion slerp are ported to C in this toolbox and verified against referenceSystem** (owner) | referenceSystem's C API offers the rotation only, and the native right-hand side cannot reach MATLAB; a copy with a mechanical comparison is R1-compliant, and upstreaming stays open |
| 2026-09-15 | **NRLMSIS 2.1 is used under NRL's academic licence, fetched and compiled at setup, never redistributed** (owner) | the same pattern as NRLMSISE-00 (O-34); the licence statement travels in THIRD-PARTY-NOTICES.md and is printed at fetch time |
| 2026-09-15 | **Accommodation: `pilinski2010` and `walker2014` are delivered, `precise` uses `walker2014`; SESAM 2013 is not implemented** (owner) | the 2013 paper is paywalled and its sticking law and shape factor are in no source opened; Pilinski's 2011 thesis prints the 2010 constant (K = 4.98e-17 m³/K, §7.3), and the adsorbate mixture is restated in open papers. The frequently repeated 7.5e-17 appears in none of the primary sources opened and is not used (R3) |
| 2026-08-27 | Concept 3 from scratch; both predecessors are input only (owner ruling 1b) | see Part C |
| 2026-08-27 | Gravity, EOM, integrator, least squares out of scope (P9) | ecosystem has owners; v0.1.1 F1 |
| 2026-08-27 | No v1 compatibility layer | F2 |
| 2026-08-27 | `sunGCRF` removed; tier is an option | F3 |
| 2026-08-27 | S0 = 1361 W/m²; 1367 only as `srp.tsi` override | F7, IAU 2015 B3 |
| 2026-08-27 | Seven test categories; exemptions reasoned | F6 |
| 2026-08-27 | Sensor layer with abstract `Sensor`; quantum and other instruments are subclasses (owner ruling 3: "also other") | P8 |
| 2026-08-27 | `qc` package as LP-5 between sensor and calibration (owner ruling: "add it") | F11 |
| 2026-08-27 | Additional effects (thruster forces, self-shadowing ray tracing, antenna thrust, Lorentz, lunar eclipse, transient thermal) are in scope as opt-ins (owner ruling 2) | never in a preset |
| 2026-08-27 | Handbook normals accepted with 1e-5 tolerance | measured: 1 − 5e-7 |
| 2026-08-27 | Speed ratio bound 1000 for the interpreted hampel loop | measured 382; interpreter/C gap, not an algorithmic defect |
| 2026-08-27 | Errors through one raise point with taxonomy and hint | F12 |
| 2026-08-27 | Run logs live in `tests/last_run.txt`, not in this document | F4 |
| 2026-08-27 | CI = audit (Python, no MATLAB) → Code Analyzer (0 findings enforced) → `ciGate`; dependencies checked out as siblings; `main` protected | referenceSystem precedent; cheapest gate first |
| 2026-08-27 | Stubs assign `varargout = cell(1, nargout)` before `fail`; sensor stub methods use `~`; `mustBeBetween` replaces `mustBeInRange` | Code Analyzer must be clean without suppression tags (which went stale: MSNU) |
| 2026-08-27 | Python mirror of the kernel layer only, frozen after LP-5 (owner ruling C1 ii); no mirror of the other layers | a full mirror lengthens the chain GROOPS↔Python↔MATLAB and inverts "MATLAB is normative" |
| 2026-08-27 | GROOPS fixtures are generated once and committed; GROOPS is never a build or test dependency (owner ruling C2 i) | reproducibility; CI stays MATLAB + Python only |
| 2026-08-27 | LP-2 before the cross-validation tooling (owner ruling C3) | both need the same real L1B day; the readers are the prerequisite of every later oracle |
| 2026-08-27 | Deliverables move by zip or git, never by base64 transcription | a transcribed bundle failed its MD5 in this session; caught before anything landed |
| 2026-09-06 | Version numbering re-based on the owner's 0.1.0 (PR #7): the numeric fields of `toolboxVersion` are the authority and `semantic` is composed from them; minor = levels delivered (LP-2 → 0.2.0), patch = HANDOVER amendment without code change; every copy (Contents.m, README, doc/html footers) is compared mechanically by `auditPackage.py`; tests compare against `toolboxVersion`, never pin a literal | R1: PR #7 left `major = 2` beside `semantic = "0.1.0"` (CLAUDE.md §6 open defect) |
| 2026-09-06 | Real-data fixtures are committed: `data/samples` holds the five 2020-01-01 GRACE-C RL04 products as the predecessor shipped them (THR1B and MAS1B the complete day; ACT1B, GNV1B, SCA1B the first 600 records with `num_records` rewritten); `data/spaceweather` holds CelesTrak SW-All and SET SOLFSMY/DTCFILE excerpts for 2019-09-01..2020-03-31. Tests never skip for absent data: a skip counts as `incomplete` and reddens the gate | no full-day 1 Hz file exists on the host (searched 2026-09-06); an excerpt of real data is real data; a skipped test is not evidence |
| 2026-09-06 | `readGNV1B` takes `'Eop'` explicitly (readEop struct or file path) and has no default series; `'Frame'` "GCRS" (default) or "ITRS" (no EOP) | R5; two EOP series differ by 9 mm on the arc, so the choice is the caller's |
| 2026-09-06 | `readSCA1B` conjugates once through `referenceSystem.invertQuaternion` and returns its canonical q0 ≥ 0 form; the sense is proven on the real sample (body +Z 0.86° from nadir; 96° when not conjugated) and re-tested every run | no hand-written quaternion algebra (A.2) |
| 2026-09-06 | `readTHR1B` returns `qc.eventTable` rows (class impulse, source THR1B, axis "all", amplitude and sigma NaN, confidence 1; records without any on-time are counted, not listed) and the per-thruster on-times in `aux`; the thruster-to-axis map is in no source in hand and is not guessed | R3; the detector (WP-31) assigns axes from the data |
| 2026-09-06 | `readMAS1B` returns gas masses only and refuses a total ≥ 100 kg; `massModel` anchors the handbook constant at the first MAS1B sample, interpolates linearly, holds the nearest value up to two sample intervals beyond the series and refuses further (`inputEpoch`); every mass it returns must exceed 500 kg | C.1 (MAS1B ≠ spacecraft mass); the first hour of a day has no MAS1B sample (records at 3599.99 s + k h) |
| 2026-09-06 | `spaceWeather` serves two families that are never mixed: celestrak (MSIS drivers: previous-day observed F10.7, centred 81-day observed mean, daily Ap, the 7-element ap history, 3-hourly Kp) and jb2008 (F10/S10 one day, M10 two days, Y10 five days earlier, 81-day values with the same lags, DTC at the UTC hour; the lags as printed in the SOLFSMY header); a query outside the table or hitting a missing value errors (`dataRange`, `dataMissing`) and never returns NaN | E.2 rule 7; a NaN driver silently zeroes or corrupts density |
| 2026-09-06 | `readMAG1B` remains a stub (O-8); `importCERES`, `igrfCoefficients` stay in LP-3/LP-3b as planned | no oracle, no unit |
| 2026-09-06 | `l1bLoad` calls `gracefoReading.readL1b` for one file and `readL1bMulti` for several (the dependency function that fits the input; the aux layouts are normalised to the plural names) | O-12: the multi-file reader crashes on an empty result; found by the `dataEmpty` tests, all six red rows of the first run were this one cause plus one test-authoring regexp |
| 2026-09-08 | **Superseded**: `l1bLoad` routes every read through `readL1bMulti`, one file or many; the branch and the fabricated `products`/`fileNames` are deleted | O-12 is fixed upstream in gracefoReading 1.0.1, so the branch's only reason is gone; one code path is R1 and makes `dataEmpty` uniformly reachable |
| 2026-09-08 | The foreign `noRecordsInRange` warning is NOT suppressed before `dataEmpty` is raised | it names the parsed record count and the active filters, which the NonGravForces error does not; foreign diagnostics propagate unwrapped here as foreign errors do, and a suppression is code that rots |
| 2026-09-08 | The `gracefoReading` minimum is declared `1.0.1`, with a PATCH component, against VERSIONING.md's MAJOR.MINOR rule | the unified path genuinely needs the patch; declaring `1.0` would state a minimum this toolbox does not run on. R3: say the measured thing, then record the deviation (O-25) |
| 2026-09-08 | Version 0.5.3, HANDOVER Revision 2.7 | patch = HANDOVER amendment without a level; rule of 2026-09-06 |
| 2026-09-07 | LP-3 scope (see Part G): msis00 via a MEX of the vendored, unmodified Brodowski C (release 20041227 + 20070727 fix, digests in `mex/nrlmsise00/README.md`), `gtd7d` default; the DOCUMENTATION §7 table is the oracle at its printed 7 digits | C3 sanctions this MEX; the table is the authors' own published output |
| 2026-09-07 | `buildMex` becomes public unit 52 (WP-05b): builds `+env/private/msisGateway.<mexext>` from `mex/`; the binary is never committed (platform-specific, silently stale); `setupNonGravForces` reports its state; `density("msis00")` refuses with `depMissing` naming `NonGravForces.buildMex` when it is absent; `ciGate` builds first | referenceSystem precedent (`.gitignore` rationale) |
| 2026-09-07 | Shadow and Earth-radiation geometry use a sphere of radius `ellipsoidParameters("WGS84").a` (6378137 m) measured from referenceSystem; oblateness and refraction are SOLAARS-CF's business (O-13) | conical model is the C.3 approximation; the sphere of the equatorial radius is the conservative (widest) shadow |
| 2026-09-07 | `eclipseFactor` returns nu in [0,1] per epoch with the lunar factor multiplied in when `'Moon'` is given; the conical formula is Montenbruck & Gill 3.4.2 with exact limits (umbra 0, annular 1 - b^2/a^2, full light 1) | C.3 |
| 2026-09-07 | `earthRadiance(backend, epoch, grid, r)` tiles the VISIBLE CAP in rings of equal central angle around the sub-satellite point and equal-azimuth sectors (`erp.grid` = [nRings nSectors]), as Knocke 1988 does; a first version with a fixed global 12x24 latitude-longitude grid missed the flux by 55 % at 490 km (measured against the pure-Python 0.25 deg integration) because the sub-satellite tile dominates 1/rho^2, and was replaced before any test was written. The zonal model needs no Earth rotation, so no EOP enters ERP; the seasonal phase is counted from 1981-12-22 00:00 UTC with period 365.25 d (the C.1 pitfall: not Jan 1); `pointSource` is the far-field Lambert sphere with a0, e0 and the Lambert phase law (measured 20-40 % off the tile sum at LEO: a fast tier, documented) | near-field weighting; CERES (Earth-fixed grids) will need the rotation and gets it in WP-27 |
| 2026-09-07 | Source list convention: `.dir` (3xMxK, body frame) points FROM the satellite TO the source; `.P` (MxK, Pa) is the radiation pressure at the satellite; `.band` `vis` for Sun and albedo, `ir` for longwave; M = 1 + 2 * tiles, constant over K, invisible tiles carry P = 0 | vectorisation over K needs a fixed M |
| 2026-09-07 | The Python oracles are pure Python (`math` only): no numpy is installed for the host's Python and installing packages is not this session's call | oracle independence needs no library; speed is irrelevant for a one-off table |
| 2026-09-07 | Version 0.3.0 for LP-3 (minor = level, P-7 as re-based); HANDOVER Revision 2.2 | rule of 2026-09-06 |
| 2026-09-07 | `eclipseFactor` evaluates the Montenbruck & Gill lens area with atan2 arcs and the separation angle by atan2(|d×s|, d·s): the textbook acos form lost 3e-11 in nu because the Earth-disc arc has its cosine within 2e-3 of 1 and the lens area is a difference of terms 1e4 times larger than itself (measured 8e-12 rotation invariance, 6.4e-10 vs the oracle; after: 2.5e-14, 5.6e-14) | the HANDOVER accuracy figure is 1e-12; found by the rotation-invariance metamorphic test, confirmed by the sharpened oracle |
| 2026-09-07 | `density` gains `'LocalSolarTime_h'`: the published NRLMSISE-00 test cases set lst = 16 h regardless of UT and longitude, so the reference table is reproducible only with an explicit override; the default remains the consistent sec/3600 + lon/15 the reference recommends | R5: the override is explicit and recorded in provenance; without it the table is unreachable and the oracle is lost |
| 2026-09-07 | The eclipse and Knocke oracles were sharpened in the same session (kink split and cosine substitution; uniform-sphere cases at 12 R and 100 R) until their own truncation was below the claim being tested; an oracle whose error dominates the comparison certifies nothing | 3.2e-7 of the first oracle version was the oracle's Simpson error at a square-root endpoint, not the model's |
| 2026-09-07 | Data audit on the owner's question (downloads allowed): SOLAARS-CF thesis fetched (open access) and transcribed; IGRF-14 coefficients fetched to `data/igrf` (SHA in its README) for WP-23; CERES EBAF/SYN1deg are behind an Earthdata login and cannot be fetched by the assistant (credentials are the owner's): O-3 stays with the owner, who drops one file into `data/ceres`; NRLMSIS 2.1, DTM2020 and JB2008 are Fortran sources (reachable) and the MinGW-w64 support package on the host ships `gfortran.exe`, so LP-3b can build them through a C gateway without a MATLAB port; PO.DAAC (full-day L1B, MAG1B) also needs a login: O-8 stays with the owner | R3 and the credential rule |
| 2026-09-07 | `solaarsCF` in `eclipseFactor`: nu from Robertson's three-sigmoid fit (eq. 3.86) of (rR, rE) in 1e6 m; rE oblateness-scaled by a/b of WGS-84 in the ecliptic-parallel component; the ecliptic normal is MEASURED from the Sun state (r_sun × v_sun) passed as `'EclipticNormal'` (radiationSources supplies it from `geocentricBodyState`'s velocity); day side (rR ≤ 0) is full light; rR clamped to the fit domain [2, 45] Mm with the clamp count reported | no obliquity constant is assumed; the fit domain is the thesis's |
| 2026-09-07 | `geom.graceFO` plate polygons are built to CLOSE on each other: the nadir half-width is where the printed 50° flank meets the nadir plane (0.9525 m, not the printed 1.9440/2), the overhang is the printed inner-flank width beyond it; polygon areas stay within 2 % of the printed areas, which remain the force authority. Centroids are the polygons' area centroids; the LP-0 centroids had the array flanks 0.39 m too high and the inner flanks 0.5 m off (found by the polygons, confirmed by a 2.5 % spurious self-shadowing of the flanks from straight above). `MacroModel` now refuses a centroid more than 5 cm from its polygon's centroid (R1) | one geometry, derived once; the handbook's own dimensions are inconsistent at the 2-4 cm level |
| 2026-09-07 | `selfShadow`: stratified ray cast on the plate polygons (default 40 x 40 samples, error ≤ ~2/m at a shadow edge) and an icosphere look-up table with barycentric interpolation (`erp`-style tier default 642 vertices); the pure-Python ray cast `tools/pyxcheck/selfshadow_check.py` (160 x 160) is the oracle; the LUT-vs-raycast difference at grazing incidence is large by construction (the table cannot resolve a shadow edge sharper than its vertex spacing) and is recorded, not hidden | T4 convergence sweep over Samples and LutResolution |
| 2026-09-07 | `AttitudeProvider` is a value class over `referenceSystem.interpolateQuaternions` (no extrapolation, foreign refusal propagated) with the static `nominal(r, v, xAlongVelocity)` (the sign is the caller's, R5: GRACE-C flew +X against the velocity on 2020-01-01); `checkAttitudeSense` reports the nadir angle for q and for q* so the verdict shows the alternative would fail | HANDOVER C.1 pitfall made mechanical |
| 2026-09-07 | `thermalState`: steady state (Wöske eq. 5.17) and the lumped transient wall (eq. 5.18) integrated by RK4 with dt ≤ τ/10; heat capacities per area from Wöske Table 6.1 stacks live on the macro model; internal power needs an explicit share (R5) | the transient tier is an opt-in; its capacities are educated guesses and say so |
| 2026-09-08 | **A test called `corr`, which is the Statistics toolbox: green on the owner's host, red on CI.** HANDOVER D.7 C2 allows base MATLAB and the PCT only, and the owner's host has Signal Processing and Statistics installed while the runner does not, so the constraint is invisible locally by construction. The correlation is now computed in base MATLAB, and `auditPackage.py` scans every `.m` file under `+NonGravForces`, `tests` and `ciGate.m` for a list of toolbox-only names (fault-injected in the same round). WP-21 had scheduled that scan for LP-R; one red CI run paid for doing it now | a constraint no instrument checks is a wish |
| 2026-09-08 | **`geom.selfShadow`'s look-up interpolation returned NaN, and no test could see it.** The face search accepted a face the ray was nearly parallel to: the ray scale explodes, the three barycentric areas come out [0 0 0], that passes `bc >= 0`, and `b/sum(b)` divides by zero. Measured at resolution 162: 8 of the first 40 icosphere vertices returned NaN for every plate. The metamorphic test read `max(abs(fv - table), [], 'all')` — which OMITS NaN — and so reported 0. Found only because a clamp added the same day turned the NaN into a 0 and the test finally failed. Fixed three ways: the face test now requires a positive denominator and weights summing to one; a non-finite fraction is refused rather than clamped away (R4); and `NgfTestCase.verifyBoundAndRecord` now fails on a non-finite measurement, so no bound anywhere in the suite can pass because a reduction dropped a NaN | an instrument's silence is evidence only about what it inspects |
| 2026-09-08 | **Count predicted 263, ran 261 on the branch, 262 after merging main.** Reconciled twice: 19 tests were added (7 wind, 7 partials, 5 Model) but `tNgfStubs` is parameterised over the stub list and shrank by 2 as `env.wind` and `kernel.partials` left it, which the prediction forgot; then main's concurrent 0.5.3 added one more. That 0.5.3 change (the O-12 discharge) had been written on Tier B and marked "NOT been run, predicted 245": this run includes its test, so its debt is discharged and its prediction was right | rule 5: a missed prediction is a finding even on green |
| 2026-09-08 | The gas force is not antiparallel to the relative velocity and the tests no longer say it is | Sentman's shear on the 50 deg flanks is genuine lift: measured 0.81 deg unshadowed, 0.34 deg with self-shadowing, against 1e-16 for constantCd, which carries no shear by construction |
| 2026-09-08 | `Options` gains `drag.cd` (the coefficient constantCd uses, 2.2) and `srp.crScale` (the radiation counterpart of `drag.cdScale`) | `Model` needs an actual coefficient for the constantCd branch, and calibration needs a radiation scale to estimate; this toolbox has no lumped reflectivity coefficient because its optics are per panel |
| 2026-09-08 | `Model` gains `'Eop'`, `'SpaceWeather'` and `'Ephemeris'` at construction | `accel` needs the geodetic height (so the ITRS rotation, so EOP) and the density drivers; a Model built without them refuses by name when evaluated, rather than assuming a file or a default driver set (R5) |
| 2026-09-08 | The relative velocity is `env.wind`'s to compute, and its co-rotation term is `v - Rdot*(R'*r)` taken from `referenceSystem` | one authority for a fact the theory document warns about (R1, theory §2.2: the co-rotation term is formed in `applyRotation`, never by hand). Measured 255.8 m/s at −59.5° latitude on the sample arc |
| 2026-09-08 | The complex-step verifier named in the LP-5 and LP-6 prompts is impossible on this chain; Richardson-extrapolated central differences replace it | measured: MATLAB's `erf`/`erfc` raise `notFullReal` on complex input and `abs`/`vecnorm` drop the imaginary part, so no complex step survives Sentman's coefficients or any norm. Recorded rather than quietly substituted (R3) |
| 2026-09-08 | LP-6b ports `thermalRecoil` too, which its prompt did not name | the frozen mirror covers it, it is six lines of C, and without it the compiled path cannot return the acceleration `Model.accel` returns — which is the only reason the path exists. Scope widened deliberately and recorded rather than taken silently (P-3). **APPROVED by the owner 2026-09-08** |
| 2026-09-08 | LP-6b does NOT link `rsApiRotation`, contradicting its own prompt and B.1. **APPROVED by the owner 2026-09-08** | `mex/rsApi.h` was read: `rsApiRotation` returns the ITRS→GCRS matrix, and the frozen kernel is body frame from end to end, so there is no rotation for it to evaluate. The prompt line survived from a planning revision written before the kernel's frame convention was settled. B.1 and D.6 are corrected; the header remains the right route for a C *integrator*, which is a caller and not this toolbox |
| 2026-09-08 | D.6's justification for LP-6b is rewritten: it is not the 5-15 us call overhead | LP-6's own O-29 measures `Model.accel` at 88 ms per epoch (standard) and 2 ms (fast) -- both superseded by WP-C3 and by the 2026-09-18 re-measurement -- where 15 us is 0.02 % and 0.75 %. A premise contradicted by the project's own measurement is not a premise (rule R3). The port keeps three other justifications, all of them stated; the speed question becomes O-31 and is answered by a number, not by an argument |
| 2026-09-09 | **Owner ruling: all data is optionally downloaded during setup; otherwise the user supplies it** | it puts data on the same footing as the NRLMSISE-00 source, which is the consistent position: the toolbox names what it needs, fetches what it can, and a reader who declines supplies it. `data/` stays in the REPOSITORY, which the ruling did not ask about — the three providers publish redistribution terms where NRLMSISE-00 granted none, and removing it would make every test and demo depend on two networks staying up. The `.mltbx` still carries none, which is what the ruling settles |
| 2026-09-09 | The documentation build is `tools/buildDocs.py`, not `tools/buildDocs.m` as the LP-R prompt named | it has to run where MATLAB does not. The audit job is the fast gate and has no MATLAB licence, and a whole session was spent on Tier B this week with the bridge down; a doc check that needs the scarce Windows runner is a doc check that runs last. Help blocks are plain text and parsing them needs no MATLAB. Folding `buildFunctionList.py` into it keeps one generator for `doc/html` rather than two (R1) |
| 2026-09-09 | **Owner ruling: track upstream NRLMSISE-00 and take its bug fixes** | the vendored 2007 copy was not reproducible from upstream (no commit in the maintainer's history carries its digest), so a public repository could not pin to it without stranding every reader. Upstream's fixes are real -- Argon densities, a memory corruption -- and its published oracle table is byte-identical, so the move costs no reference value. CI pins a COMMIT rather than a branch: the ruling is to track upstream deliberately, not to let someone else's push change this toolbox's densities silently |
| 2026-09-09 | **Owner ruling: MIT for the toolbox's own code** | recorded in `LICENSE`. The licence text itself carves out `mex/nrlmsise00/`, because a blanket MIT over vendored source this project has no grant for would be a claim it cannot make. The vendored question becomes O-34 and `THIRD-PARTY-NOTICES.md` |
| 2026-09-09 | **Owner ruling: O-29 accepted — 88 ms per epoch (a superseded figure) is not a release blocker** | it is recorded as a known cost of the standard preset rather than a gate on LP-R. The fast preset is 2 ms, and O-31's measurement already says the fix would be a `geom.selfShadow` pass and not the C kernel, so shipping now costs nothing that a later profiling session cannot recover |
| 2026-09-18 | **Owner ruling: the `tools/` smoke test stays, and its gate time is accepted** | `tests/tNgfTools.m` runs `tools/benchmark.m` under the gate and costs about 10 s (192.3 s against 182.7 s), almost all of it reading the EOP file and building the two Models rather than timing anything. Asked because the instrument is explicitly NOT gated by design -- it measures the machine as much as the code -- and the answer is that RUNNING it is not the same as GATING its numbers: completion is asserted, no duration is. Same trade the project already accepts for `tNgfDemos` |
| 2026-09-18 | **Owner ruling: O-6 is parked, "open until further notice"** | it stays on the books as a recorded gap and stops being an ask. No session re-raises it, schedules it, or counts it against LP-R until the owner lifts the parking. The standing request -- one GROOPS `SimulateAccelerometer` run for GRACE-C 2020-01-01 with per-effect output, which would close O-19 and WP-34 with it -- is recorded in the row and is not to be repeated as a question |
| 2026-09-09 | **Owner ruling: lsqvce is overkill for WP-19; `app.calibrate` carries the toolbox's own solver** | the problem is linear and separable by axis — three QRs of a few dozen unknowns. Variance components, spectral models, streaming accumulation and arc factors have no work to do on it. The dependency row is REMOVED rather than kept as optional-but-unused, because a table naming a toolbox nothing calls will be believed. O-24 and O-33 are void with it, and the CI runner needs no new checkout or token scope. What was given up is honest formal sigmas under coloured noise, which `calibrate` states in its own help |
| 2026-09-08 | `sensor.Electrostatic` fixes its own omega sign convention in its help instead of deferring to `referenceSystem` | the two do not obviously agree: that unit's R runs the other way and its sense is configurable, and the disagreement is invisible in the centrifugal term. A measurement model whose sign depends on an unread convention is a model that will be wrong by exactly -1 in one term and pass most tests. The contract is stated as `skew(omega) = R' Rdot`, derived here, and the relation to referenceSystem is recorded from a run rather than guessed (O-32) |
| 2026-09-08 | The CRN filter is an argument, not an implementation | Wu et al. (2006) defines the kernel and has not been read for this work. A guessed low-pass would make model and L1B comparable in appearance and not in bandwidth, which is worse than no filter because it looks right. The class applies any FIR kernel, marks the (L-1)/2 edge samples invalid rather than faking them, and filters the partials identically |
| 2026-09-08 | The gravity-gradient term is an optional caller input, not a computed one | T = d(g)/d(r) needs a gravity field and this toolbox has none by design (D.7 C4, R1). Omitting it silently would be an approximation the class makes on the caller's behalf; requiring it would make the class unusable without a gravity model for a term of order 1e-12 m/s^2. It enters where the observation equation puts it when supplied, and `y.info.gradient` says which happened |
| 2026-09-08 | **P-11 (owner ruling): the vectorised MATLAB path and the pointwise C path are BOTH first class. Neither replaces the other and neither is justified by being faster than the other.** | `Model.accel`/`series` vectorise over K epochs and are the path for an arc; `kernel.singleEpoch` evaluates one epoch in C and is the path for a right-hand side called once per step. They answer different questions, so the speed ratio between them is a measurement and not a decision (O-31 drops to informational). What the ruling DOES oblige: two paths that are both kept must be kept in step, which is `tests/tNgfMex/singleEpochEndToEndPrecision` — the compiled path composed over the real arc must return the acceleration `Model.accel` returns for that epoch |
| 2026-09-08 | `Model` is not routed through the compiled kernel | under P-11 there is nothing to route: `Model` IS the vectorised path. A `'Backend'` option on `Model.accel` would make one path a mode of the other and put a fidelity-tier fallback where D.2 forbids one. A caller that wants the pointwise path calls `kernel.singleEpoch` |
| 2026-09-08 | `buildMex` returns a 1xT struct array, one element per MEX target | a second gateway had to go somewhere, and a scalar report can only describe one. `'Target'` names one and keeps the scalar shape, which is what the env tests and the robustness test now pass. `ciGate` and `setupNonGravForces` loop. The alternative — a second function `buildKernelMex` — would have duplicated the compiler search, the digest scheme and the staleness rule, i.e. three second authorities (R1) |
| 2026-09-08 | The C is gated by CI in the `audit` job, not the `matlab` one | `ubuntu-latest` has a compiler and the self-hosted Windows runner is the only machine in the project with a MEX one, so a C typo would otherwise wait for the expensive job on a scarce machine. `mex/test/runCChecks.sh` compiles the gateway against a `mex.h` shim (no link) and runs the mirror parity harness; both need neither MATLAB nor numpy, which is the same constraint every oracle generator here already meets |
| 2026-09-08 | `app.forceBudget` (WP-30) moves from LP-6 to LP-7 | it is the application layer, where `calibrate` and `densityInversion` already live, and it is a thin reduction over `Model.budget` which lands in LP-6. `Model.m`'s own stub had said a later level all along, so Part G and the code disagreed |
| 2026-09-07 | **Count predicted 243, ran 244.** The miss is mine: six tests were added across four classes over the course of the sweep and the prediction was written after five of them. The reconciliation is recorded because a missed prediction is a finding even on green (rule 5): tNgfCore +2, tNgfQc +2, tNgfGeometry +1, tNgfEnv +1. Separately, the file-side counting script used until now over-counted, because the regexp `methods \(Test[^\n]*\)` also matches `methods (TestClassSetup)`; MATLAB's own suite is the authority | rule 5 |
| 2026-09-07 | `qc.eventTable` accepts a NaN amplitude and sigma (meaning NOT ESTIMATED) and refuses only an Inf, while confidence must be finite | found by running: the first, blanket refusal broke all seven `readTHR1B` tests and one `detectImpulses` test, because a thruster firing is a certain epoch of unknown size and `readTHR1B` writes NaN there deliberately. The real defect was only ever the probability |
| 2026-09-07 | `Model` validates its `options` argument itself instead of through the `arguments` block | the block raised `MATLAB:validation:UnableToConvert`, while the help had documented `NonGravForces:Model:inputOptions`, which nothing could raise; the test asserted the undocumented identifier and was corrected with the code |
| 2026-09-07 | The shortfall filter of the gate is `~complete`, not `~complete & ~exempt` | measured: the second term hid three real coverage holes (`const` and `toolboxVersion` had no robustness test, `eventTable` had neither robustness nor metamorphic) behind exemptions granted for OTHER categories, and every run since LP-0 had printed 0 audit shortfalls. `auditTestCategories` already credits exempted categories into `complete`, so the term was never needed. The four tests were written rather than the exemptions widened |
| 2026-09-07 | `qc.eventTable` refuses non-finite amplitude, sigma and confidence | found by writing the missing robustness test: a NaN confidence passed `any(cf < 0 | cf > 1)` silently, because every comparison with NaN is false (R4) |
| 2026-09-07 | `auditPackage.py` gains two checks, each fault-injected in the same round: the catalogue must read as seven comma-separated fields, and the vendored NRLMSISE-00 sources must match the digests their README records | an unquoted comma in `env.eclipseFactor`'s h1 had truncated its purpose to "Shadow function nu in [0" in Contents.m and FunctionList.html, and the vendored digests were recorded but compared by nothing. Both checks were verified to fire on a broken tree and stay silent on a healthy one |
| 2026-09-07 | `runNgfTests` refuses a `'Name'` pattern that selects no class instead of silently running the whole folder | a mistyped pattern would otherwise report a green run of tests the caller did not ask for |
| 2026-09-07 | The obligations digest in CLAUDE.md is replaced by a pointer to Part H | it had drifted to O-15 of 21 and showed two closed items as open, in the one file the harness injects into every agent; a duplicate that no check compares is worse than a reference (R1) |
| 2026-09-07 | The independence of the mirror is claimed from ARITHMETIC, not from access control: `CLAUDE.md` reaches every subagent whatever the prompt says, and the self-reported contamination flag was inconsistent (2 of 14 set it while naming only permitted sources) | measured during the run; recording the leak is worth more than a claim of purity that the harness cannot honour |
| 2026-09-07 | The CLAUDE.md pitfall "Sentman: the shear term divides by s^2, not s" is corrected to name its coefficient and convention: in this toolbox's normal/tangential convention the TANGENTIAL coefficient is l[E/(s sqrt(pi)) + gamma Z] and contains no 1/s^2 term; the lone Z/(2 s^2) belongs to the LIFT coefficient of the drag/lift pair | three independent derivations contradicted the note and the edge-on limit settles it (ctau = l/(s sqrt(pi)) exactly at theta = pi/2, measured to 5e-16 on both sides). The Part C D1 audit row is a historical claim about the predecessor's code and is annotated, not rewritten |
| 2026-09-07 | The mean-molar-mass speed ratio of the reference document (section 2.3) is an approximation, not an equivalent of per-species mixing; both drag derivations rejected it independently and measured the error at the 490 km composition: harmless on a ram face, -24 % at 100 deg incidence, total collapse at 135 deg | the kernel already mixes per species; the document is to be corrected at LP-R with WP-21 |
| 2026-09-07 | `hampel_check.py` recomputes the committed flags from the committed series instead of redrawing it | the series is data (like the L1B excerpts); what must be re-runnable is the algorithm. It had imported numpy, which no tier has, so the one oracle in the tree could not be regenerated at all |
| 2026-09-07 | The completeness-critic stage of the mirror workflow died on a session limit; the audit was performed by hand and passed: stdlib-only imports in all six modules, all six import together, and boundary cases run (theta = 0, pi/2 and pi; alpha = 0 and 1; empty plate lists; zero density; an exactly edge-on plate; a source exactly in-plane; T = 0 K; r = 0, which refuses with a reason) | a stage that did not run must be named, not silently dropped |
| 2026-09-07 | LP-5x takes a PATCH bump (0.5.1), not a minor one | the version rule maps minor to the level delivered, but LP-5x adds no public unit and changes no interface: it is cross-validation tooling plus a HANDOVER amendment. 0.6.0 stays with LP-6, where the Model layer lands |
| 2026-09-07 | WP-33 is a standard-library mirror, not a numpy one | measured: the host has no numpy and the ubuntu audit job installs nothing, so a numpy mirror is regenerable on neither tier; `hampel_check.py` had already become unrunnable for exactly this reason |
| 2026-09-07 | The mirror's deriving agents were denied the kernel sources, the existing Sentman oracle, the kernel tests, HANDOVER.md and CLAUDE.md; the plate tables are given to them as data | T9: an oracle transcribed from the implementation inherits its misunderstandings. The records were denied too because the LP-5 decision rows quote implementation detail (erfc, the accommodation bracket) |
| 2026-09-07 | Version 0.5.0 for LP-5; HANDOVER Revision 2.4 | rule of 2026-09-06 |
| 2026-09-07 | Version 0.4.0 for LP-4; HANDOVER Revision 2.3 | rule of 2026-09-06 |
| 2026-09-07 | LP-5 kernel returns normal and tangential plate coefficients (Cp, Ctau) derived from the momentum integrals, not the drag/lift pair of Doornbos (2011) that the reference document mislabels; the tests assert both projections | the hyperthermal limit exposes the mislabel (2 gamma vs 2 gamma^2); the numerical oracle is independent of either closed form |
| 2026-09-07 | `accommodation` delivers fixed only; sesam/langmuir refuse `notImplemented` (O-16); `gsiCoefficients` refuses cll/schamberg with `unsupportedModel` (O-17) | R3 measured never declared: no constant is transcribed from memory |
| 2026-09-07 | `qc.residuals` takes the preliminary calibration as explicit name-values ('Bias', 'Scale') instead of a Sensor object; the Sensor route is added with LP-7 through the same unit | sensors are LP-7; a placeholder interface would be a second authority |
| 2026-09-07 | Upstream defects filed as issues: referenceSystem#15 (O-10), gracefoReading#6 (O-12) | owner request of 2026-09-07 |
| 2026-09-07 | O-2: `MacroModel.N` holds renormalised unit normals, `panels(k).normal` the printed digits; kernels renormalise nothing | one renormalisation at the single validation point (R1); the printed six digits are a truncation of a unit vector (norm 1 − 7.3e-7, measured), 7e-7 in every per-plate force otherwise |
| 2026-09-07 | `detectImpulses` correlates with the zero-mean template (h − mean h) | a 4-sigma level step produced a plateau of false pulses with the raw template (measured on the synthetic case); the zero-mean filter is blind to offsets and steps except at their edge |
| 2026-09-07 | `detectJumps` BIC penalty is (p + 1) sigma^2 log K (segment parameters plus the change-point location, the changepoint-package convention), not p sigma^2 log K | measured: 1000 white samples with one 4-sigma step gave 4 change points under p log K, 1 under (p + 1) log K; the exhaustive oracle is penalty-agnostic, the penalty is a statistical claim |
| 2026-09-07 | `gsiCoefficients` evaluates Z = erfc(-gamma s) and the re-emission bracket as (1 - alpha) + alpha 2Tw/(Tinf s^2) | 1 + erf(gamma s) cancelled to 0 or below on leeward plates (negative Cp of -9e-20 measured at s = 8, theta = 143 deg); the bracket 1 + alpha(x - 1) lost 1e-8 relative at alpha = 1 (measured 4e-13 absolute against the closed-form bound) |
| 2026-09-07 | `thermalRecoil` accepts T = 0 K (emits nothing) | `thermalState` steady state puts unlit plates without internal power at exactly 0 K; refusing them broke the steady-state metamorphic test (found by running) |
| 2026-09-07 | `qc.screen` replaces only ISOLATED outliers (runs of <= 2 flagged samples) before the matched filter, all outliers before the jump search, and lets detected pulses claim their samples | Hampel flags the peaks of a 6-sigma pulse in a run of ~6 samples; median-replacing them turned the pulse into four level jumps (measured); a spike is 1-2 samples |
| 2026-09-07 | PELT in `detectJumps` keeps not-yet-admissible candidates (segments shorter than MinSegment) in the candidate set instead of pruning them by omission | first implementation returned no change point at all on a 5-sigma step: the candidates were dropped before they became admissible; found by the exhaustive oracle, not by reading |
| 2026-09-07 | `checkAttitudeSense` and the tests take angles as atan2(|a×b|, a·b): acos of a dot product near 1 loses sqrt(eps) (1.2e-6° measured on an exact nadir frame; 2e-14° after) — the same lesson as the eclipse lens area | conditioning, second occurrence in two sessions: make it the house rule for angles |
| 2026-09-07 | `selfShadow` sums the sample weights column-wise (not through a BLAS product whose summation order depends on the batch size) and sets fully covered / fully exposed plates exactly: a vectorisation test is bitwise by design | T1 vectorization category |
| 2026-09-06 | The local `gracefoReading` clone was found at a broken initial commit (`parseDataSection.m` with a truncated signature) and fast-forwarded to its origin `main` (v1.0.0) before any measurement | dependencies are used at their published state, never patched locally |

---

## Part J — The nine-expert review (positions, agreements, simplifications)

The design was argued by nine roles; the outcome is Parts C–F. Recorded here are the points where the roles disagreed and how it was settled.

1. **Thermosphere physicist** (density, GSI): wanted CLL and SESAM in `standard`. *Settled:* opt-ins; the 10–30 % density error dwarfs the GSI refinement; `precise` gets SESAM.
2. **Radiation-pressure specialist**: wanted TRP transient in `precise` and SOLAARS-CF everywhere. *Settled:* steady-state TRP in `standard`/`precise`, transient opt-in; SOLAARS-CF in `precise` only (conical is C⁰ at the umbra edge but adequate for drag work; the C¹ requirement is documented for integration).
3. **Accelerometer/calibration expert**: insisted the CoM-coupling term is a sensor property, not physics; insisted on jump columns rather than bias-spline absorption. *Adopted* (P8, `qc`).
4. **Geodesist (frames/time)**: refused any second authority for time/frames/Sun. *Adopted* (F3, F8).
5. **Quantum-sensor physicist**: the CAI needs true acceleration plus rotation environment, dead-time-aware sampling and an ambiguity model; hybridisation with the electrostatic channel is the realistic operational mode. *Adopted* as `AtomInterferometer` with an `omega` input and a `Hybrid` composition planned in LP-7.
6. **Software architect**: one raise point, one schema per artefact, value classes, no wrappers, no compat layer, one dispatcher per model family. *Adopted* (P10, F13).
7. **MATLAB performance engineer**: K-vectorised 3×K layout, `parfor` only over chunks, no GPU, MEX only where the interpreter is the bottleneck (density inner loop, single-epoch RHS). *Adopted* (D.6); measured the hampel loop gap (O-1).
8. **Test engineer**: metamorphic relations before oracles, independent Python oracles, mechanical category audit, ratios not absolutes. *Adopted* (T8, T9); caught the two authoring errors of this session by exactly this discipline.
9. **Data-processing expert (L1B)**: THR1B/MAS1B/eclipse as *priors* for the detectors, never re-discovered; ACT vs ACH product awareness; CRN response as the matched filter. *Adopted* (LP-5).

Simplifications accepted by all nine: 51 units instead of v0.1.1's ~70 planned; three presets; one kernel path for all radiation; sensors replace instrument-specific calibration code; the toolbox is smaller than either predecessor and complete.
