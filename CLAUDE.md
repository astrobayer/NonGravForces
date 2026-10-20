# CLAUDE.md — working notes for NonGravForces

Orientation file for an AI assistant working in this repository.
`HANDOVER.md` is the single source of truth; this file is the short form.
Read `HANDOVER.md` Part A, then the level prompt in Part G for the level being worked.

Created 2026-09-06 by Claude Opus 5 (Anthropic); updated 2026-10-18 (release 0.16.0) by Claude Opus 5.5. Toolbox version 0.16.0.

---

## 1. Purpose

MATLAB toolbox that computes **non-gravitational forces and torques** on GRACE-FO
(and, by construction, other LEO satellites): atmospheric drag, solar radiation
pressure, Earth radiation pressure (albedo + longwave), thermal recoil, and the
corresponding torques. Primary use case: accelerometer data analysis and
calibration for GRACE-C/D.

Successor to two rejected predecessors (`SatelliteNonGravToolbox` v1,
`NonGravToolbox` v0.1.1); redesigned from scratch after an audit found 11
quantitative defects in the first and 14 architectural flaws in the second
(`HANDOVER.md` Part C).

**The rule that governs every other rule** (reworded 2026-10-15 under N-06; HANDOVER A.2):
> NonGravForces computes forces and torques, and applies them; it implements nothing it can call.

No gravity field, no equations of motion, no integrator, no least-squares code. The `app` layer
composes (calibration through `timeSeriesPreparation.run` and `lsqvce.fit`); it does not solve.

---

## 2. Architecture

Eight layers; a layer calls only lower layers plus the four dependencies. L7 is vacant: quality
control moved to `timeSeriesPreparation.qc` (N-01) and its forwarders were deleted in 0.14.0 (O-47).

```
L8 app      calibrate, fitCalibration, compareCalibrations, perRevolution, attributeRevolutions, plotCalibration, plotAttribution, Calibration, syntheticPair, forceBudget, densityInversion, DynamicsSource (gjToolbox entry, native C)
L6 sensor   Sensor (abstract), Electrostatic, AtomInterferometer, DragFree
L5 Model    accel | torque | series | partials | budget   <- single point of access
L4 kernel   gsiCoefficients, accommodation, dragPanelForce, radiationPanelForce,
            thermalRecoil, assemble, bodyTorques, partials, singleEpoch (C/MEX)
L3 geom     MacroModel, graceFO, grace, boxWing, selfShadow, flowShadow, thermalState,
            AttitudeProvider, checkAttitudeSense
L2 env      density, wind, eclipseFactor, earthRadiance, magneticField, radiationSources
L1 io       readACC1B/SCA1B/GNV1B/THR1B/MAS1B/MAG1B, massModel, spaceWeather,
            importCERES, importGroopsGrid, igrfCoefficients
L0 core     const, toolboxVersion, dependencyReport, fail, provenance, Options, buildMex
```

59 public units, catalogued in `doc/functionlist.csv` (the authority; CI checks
both directions against the files on disk, and this sentence's number against it
— which it genuinely did not until 2026-09-16, when the sentence said 53 and the
catalogue held 55; `auditPackage.py` check 5 now compares both this number and
this file's version stamp).

### Dependencies — called by qualified name, never copied, wrapped or aliased

(lsqvce re-entered with N-02b on 2026-10-15: `app.fitCalibration` solves through `lsqvce.fit`. N-06's
records — A.2 rule, B.4 — follow with N-02c.)

| Toolbox | Owns |
|---|---|
| `referenceSystem` ≥ 2.3 (2.3.0; raised from 2.1 by N-02c for `epochTo(..., "gpsNanoseconds")`) | time scales, EOP, ITRS↔GCRS with derivatives, geodesy, quaternions, attitude interpolation, Sun/Moon states |
| `gracefoReading` ≥ 1.0.1 (1.0.1) | L1B ASCII reading, quality flags, GROOPS instrument export/import. The PATCH minimum is deliberate and is the one exception to VERSIONING.md's MAJOR.MINOR rule (HANDOVER O-12, O-25) |
| `timeSeriesPreparation` ≥ 1.4 (1.4.0) | The L1B CRN kernel (`crnKernel`), and after N-01 the `qc` stage. **Optional for the toolbox, required by `ciGate`**: `sensor.Electrostatic` filters with `crnKernel("ACC")` by DEFAULT (N-04, owner ruling 2026-10-14) and raises `Electrostatic:depMissing` without it. Its version is the `% Version` line of its `Contents.m` |
| `lsqvce` ≥ 6.35 (6.35.0; raised from 6.22 on 2026-10-19 for `R.contributions` / `solveContributions`) | The one least-squares estimator (CONCEPT Q2). **Optional for the toolbox, required by `ciGate`**: `app.fitCalibration` declares kind A (one system, one day, one arc per contiguous run) per axis its day bootstrap solves the fit's per-day contributions and its single-arc bootstrap refits kind B. Version from `lsqvce.util.toolboxVersion` |
| `gjToolbox` (odeSuite) ≥ 0.1 (0.1.0) | Gauss-Jackson integration. **Optional for the toolbox, required by `ciGate`**: `app.DynamicsSource` subclasses `gjToolbox.DynamicsSource` (so it cannot even load without it) and the `buildMex` target `rhs` compiles against `csrc/gj_native.h` (ABI 1). `Model.rhs` is gone — O-4 closed by removal, 2026-09-16 |

### Key design decisions

- **P2 one kernel**: per-panel forces are computed once; force and torque are
  reductions over the same panel loop.
- **P3 sources, not effects**: `env` produces a source list (Sun, Earth tiles,
  atmosphere); the kernel does not know which effect it is computing.
- **P4 body frame inside**: everything internal is body frame + radians + SI.
- **P5 configure once, evaluate many**: `Options` and `Model` are built once,
  then evaluated over epochs; vectorised over epochs throughout.
- **P8 sensors are adapters**: the physics core returns true dynamics
  (a_ng, τ, ω, ω̇); `sensor.*` maps them to an observable (CoM coupling, CRN
  filter, frame). Quantum (atom-interferometer) and drag-free sensors are
  first-class, not retrofits.
- **P10 one raise point**: `NonGravForces.fail(unit, mnemonic, message)`.
- **P11 two evaluation paths, both first class** (owner ruling 2026-09-08): the
  vectorised MATLAB path (`Model.accel`, `series`) for an arc, and the pointwise
  compiled path (`kernel.singleEpoch`) for a right-hand side. Neither replaces the
  other, neither is justified by speed, and `Model` gains no `'Backend'` option —
  that would make one a mode of the other. The price is that they must agree, and
  `tNgfMex/singleEpochEndToEndPrecision` is where that is checked. **LP-C adds a
  third, at L8 and not in `Model`**: `app.DynamicsSource`'s native entry
  (`mex/ngfRhs.c`) evaluates the WHOLE chain in C over tables built once per run,
  for an integrator's inner loop. Same rule — MATLAB is normative, parity is
  asserted in six configurations in `tNgfIntegration`, and `Model` is untouched.
- **MATLAB is normative.** The C/MEX port (LP-6b, hand-written, no MATLAB Coder,
  `mex/ngfKernel.c` + `kernel.singleEpoch`) and the standard-library kernel mirror
  (WP-33) are verified *against* MATLAB, never the reverse. The C is written in the
  same association order as the MATLAB lines it ports, on purpose: that is why the
  two agree to rounding instead of to "about 1e-12".

---

## 3. MATLAB conventions (agreed, enforced)

- **Errors**: identifier `NonGravForces:<unit>:<mnemonic>`, built only by
  `NonGravForces.fail`. Mnemonic taxonomy (first token): `input*`, `dep*`,
  `data*`, `numeric*`, `unsupported*`, `notImplemented`. Every message ends
  with ` -> <what the caller should change>`. Foreign errors
  (`referenceSystem:*`, `gracefoReading:*`) propagate unwrapped, never re-badged.
- **Contracts**: `arguments` blocks on every public unit. **`inputParser` is
  banned** (CI checks). Semantic checks that need several inputs go through
  `fail` right after. Nothing is validated twice.
- **Never degrade silently**: no fallback between fidelity tiers; unknown option
  fields error; a missing optional dependency errors at the point of use.
- **Style**: `camelCase` functions, `PascalCase` classes, one public unit per
  file, SI with unit-suffixed names at interfaces (`_m`, `_K`, `_W`), radians
  internally, 3×K column layout, no globals, **no `addpath` in package code, no
  printing from library code** (only `dependencyReport` may print).
- **Help template** (mandatory, CI-enforced): H1 · syntax · description · INPUT
  ARGUMENTS (Name Size Type Units Description) · NAME-VALUE ARGUMENTS with
  defaults · OUTPUT ARGUMENTS · ACCURACY (every approximation named and
  quantified) · ERRORS (identifiers) · REFERENCES with equation numbers ·
  EXAMPLE · See also · KEYWORDS · `Created <date> by <model>`.
- **Tests**: seven categories per public unit, counted by
  `auditTestCategories`: `contract`, `reference`, `precision`, `vectorization`,
  `robustness`, `metamorphic`, `speed`. Exemptions live in
  `testCategoryExemptions.m` with one reason each; `contract` is never exempt.
  A passing bound leaves its measured number in `tests/last_run.txt`.
  Assert ratios, not host-dependent absolutes. Tests assert identifiers, not
  messages. **An oracle is independent code** (`tools/pyxcheck/*.py`, published
  tables, real data), never a copy of the unit under test.
- **Only base MATLAB + Parallel Computing Toolbox.** No other toolboxes on
  default paths. No MATLAB Coder.
- **Nothing is marked ☑ without a read MATLAB run report.** `◐` =
  delivered-untested, `☑` = verified. Predict the test count before a run and
  reconcile after. **A failing test is decomposed before any bound moves.**

---

## 4. How to run

```matlab
run('...\referenceSystem\setupReferenceSystem.m')
run('...\gracefoReading\setupGraceFOReading.m')     % must run, or reader tests error
run('...\odeSuite\gjToolboxSetup.m')                % must run, or tNgfIntegration errors
cd('...\NonGravForces'); setupNonGravForces          % adds root, tests/, tools/; registers doc
runNgfTests                                          % 'ExcludeSlow',false adds the speed tag
doc NonGravForces
```

`bash mex/test/runCChecks.sh` is the C-side gate and needs no MATLAB: it compiles
the MEX gateway against a `mex.h` shim and runs `mex/ngfKernel.c` against the frozen
mirror fixtures. CI runs it in the `audit` job.

`matlab -batch ciGate` is the same gate CI runs (installs dependencies from
sibling folders, runs the suite, judges, raises `NonGravForces:ciGate:red`).

**CI** (`.github/workflows/ci.yml`): job `audit` on `ubuntu-latest`
(`auditPackage.py`, stdlib only: catalogue ↔ files, help template, `fail`
taxonomy, banned constructs) → job `matlab` on the **self-hosted Windows
runner** (checks out `astrobayer/referenceSystem` and `astrobayer/gracefoReading`
as siblings, MATLAB probe from `astrobayer/ci-actions/matlab-probe`, `codeIssues`
must be empty, then `ciGate`). Every third-party action is pinned to a commit SHA.
`main` is protected by a **ruleset** (not classic branch protection — the classic
`/branches/main/protection` endpoint answers "Branch not protected" and reading only that
would say `main` is wide open): `required_status_checks` on both jobs with strict, plus
`deletion` and `non_fast_forward`, no bypass actors. **A pull request is NOT required by any
rule** — this file claimed it was until 2026-09-08. Changes arriving by PR is a convention
here, not a guardrail. A second ruleset protects `~ALL` tags against deletion and
non-fast-forward, which is what makes a tag genuinely immutable.
Work goes on a branch → PR → green → merge → `git pull` on the host.

---

## 5. Current status (2026-10-19, day bootstrap on lsqvce contributions)

- **2026-10-19, lsqvce#183 used:** `fitCalibration`'s day bootstrap is `lsqvce.solveContributions` on the drawn
  multiplicities (lsqvce >= 6.35). 10 days of GRACE-C, 100 replicates: `calibrate` 1410 -> ~340 s, sigmas identical
  to ten digits; under 'Noise' estimate / 'Robust' the weights are now fixed per replicate (Part I, ratified 2026-10-20).
  ☑ on Tier A: 437/437, predicted 437.
- **2026-10-18, release 0.16.0:** N-05, N-09, N-11a, N-10a, #96 and the flow-shadow table released. Gate 435/435 on
  Tier A, predicted 435; benchmark re-measured (standard 0.440 ms at K=100, was 0.836); the month baseline
  passed at 0.16.0, worst 0.013 sigma, in 3.1 h (was 4.0). Issues #59/#60/#62/#66 closed; O-35 closed.
- **2026-10-17, flow-shadow table:** `geom.flowShadow` Method "table" (trilinear over the stream direction about
  +-x and s_O, exact fallback outside, flagged) is the standard preset's default (`geom.flowShadow`), read by
  `Model` and by `mex/ngfRhs.c` alike; drag moves <= 6.7e-13 m/s^2 on ten real days, ~100x cheaper per epoch.
  A 10-day `calibrate` (C): 650 -> 350 s at NBoot 0, 1541 -> 1204 s at NBoot 100 (the bootstrap is lsqvce#183).
  Tier A: 435 passed, 0 failed, 0 incomplete, GATE TRUE, predicted 435.
- **2026-10-16, occlusion speed-up:** `selfShadow`'s look-up, `radiationPanelForce` and `flowShadow`'s species
  loop restructured with every value bitwise unchanged (accel, info and torque compared with `main`):
  `Model.accel` 1.004 -> 0.85 ms per epoch, interleaved A/B. The benchmark table was not re-committed (the
  host was loaded; untouched paths moved too). Further gains change numbers: HANDOVER lists three options.
  ☑ on Tier A: 427/427, predicted 427.
- **2026-10-16, N-10a (#67):** `tools/regressionBaseline.m` calibrates 2020-01 C and D from the DATAPOOL
  (`standard`, joint over days, 100-replicate day bootstrap) and freezes `tools/baseline_reference.json`
  (C: S_xx 0.7620 +- 0.0084; D: 0.7694 +- 0.0097, 19..21 Jan without records). ~4 h, run at every RELEASE
  on Tier A: audit check 14 refuses a version the baseline was not verified at. The month is in
  `E:\DATAPOOL\GRACEFO\L1B\RL04\2020`; CI sees the pool through `READING_DATA_POOL` (#92).
- **2026-10-15, N-11a (#68, first of two PRs):** 'Regressors' in `fitCalibration` (measured series, linear onto
  the raw grid, demeaned, filtered like the bias; common over days) carried by `Calibration` (schema "3";
  `apply` refuses without the series). Measured: GRACE-FO TSU temperatures span 0.022 C per day (thermal
  control), TICUR is dead at -273.2 C. ☑ on Tier A: 427/427, predicted 427. Next: N-11b (temperature-
  dependent full scale).
- **2026-10-15, N-09 (#66):** `fitCalibration` and `calibrate` take Arc arrays: one COMMON scale, per-day biases
  as lsqvce 'day' locals (pre-eliminated), uncertainty by a bootstrap over days; 'ScaleMode' "perArc" = N
  single fits. Real fixtures: 30 min of ACT1B/GNV1B/SCA1B for 2020-01-01..03 C in `data/samples/multiday`.
  ☑ on Tier A: 421/421, predicted 421. Unblocks N-10 (the month baseline; needs a runner `DATAPOOL`).
- **2026-10-15, N-05 (#62):** `app.plotCalibration` (scale/bias with bands, residual ACF before and after AR
  whitening) and `app.plotAttribution` (per-revolution scales vs time and drivers, verdict banner); tests
  assert the figures draw their source bitwise. Found: `attributeRevolutions` needs >= 20 revolutions
  (documented; 'BlockLength' added). 59 units. ☑ on Tier A: 413/413, predicted 413.
- **2026-10-15, release 0.15.0 (O-49 closed):** N-03 released; benchmark re-run on Tier A (standard 17.0 ms at
  K=1, 0.836 at K=100), 405/405 on Tier A, predicted 405. Also on `main`: #86, a failed kernel download now
  names its cause in ciGate's message. Open N-issues: N-05 (figures, needs N-03 -- now unblocked), N-08..N-12.
- **2026-10-15, N-03c (#60, third of three PRs):** `app.perRevolution` (bias + diagonal scale per window through
  `fitCalibration`; a gated window is NaN with its reason) and `app.attributeRevolutions` (drivers beta',
  eclipse fraction, F10.7/Ap, effect amplitudes, time -> `timeSeriesPreparation.regressDrivers`; verdict at
  3 bootstrap sigma). 57 units. ☑ on Tier A: 405/405, predicted 405. **N-03 is complete; release 0.15.0 is
  O-49.** Every decision the owner had open was ratified on 2026-10-15 (HANDOVER Part I).
- **2026-10-15, N-03b (#60, second of three PRs):** `app.compareCalibrations` ranks candidate models per axis by
  blocked CV (`timeSeriesPreparation.crossval`), AIC/BIC as IC_unreliable columns; `fitCalibration` gains
  `'Span'` and `info.prediction`. Fixed an N-03a defect (constant bias + held S_yy). 55 units. ☑ on Tier A:
  390/390, predicted 390. Next: N-03c (per-revolution and attribution), then the release.
- **2026-10-15, N-03a (#60, first of three PRs):** time-varying scale. `app.fitCalibration` takes
  `'ScaleBasis'`/`'ScaleKnotSpacing_s'`; `app.Calibration` carries `scaleCurve` with `scaleAt`, `scaleSigmaAt`,
  `drift` (Wald test) and `isDriftSignificant`, schema "2"; `syntheticPair` plants `'ScaleDrift'`. A knotted
  scale finer than the bias is refused; `drift` WARNS on the formal white covariance (it flagged a planted-zero
  drift at p = 1e-18; the bootstrap gave 0.72). ☑ on Tier A: 383/383, predicted 383. Next: N-03b (model
  comparison), N-03c (per-revolution and attribution).
- **2026-10-15, release 0.14.0 (O-47 closed):** the six deprecated `qc` forwarders are deleted with
  `tests/tNgfQc*.m` (39 tests), their catalogue rows, hushes and exemptions, and the hampel/PELT oracles
  (timeSeriesPreparation carries its own); L7 is left vacant; 54 units; quality control is
  `timeSeriesPreparation.qc`. MINOR at 0.y.z per VERSIONING.md. ☑ on Tier A:
  375/375, predicted 375 (414 minus MATLAB's count of the two deleted classes, 11 + 28); benchmark re-run
  on a loaded host (compiled path 0.497 ms vs 0.457).
- **2026-10-15, release 0.13.0 (O-48 closed):** version bumped; `tools/benchmark.m` re-run on Tier A and
  `tools/benchmark_reference.csv` committed; HANDOVER's title (it said 0.11.0) now compared by audit check 5.
  Fixed cost per call rose ~2.5 ms in both presets, marginal cost unchanged (O-48 row). **The host's own
  `referenceSystem` gateway is stale** (built 2026-08-05; setup's rebuild fails while a MATLAB process
  maps it): restart MATLAB, then `setupReferenceSystem`. Next: O-47 (delete the qc forwarders).
- **2026-10-15, N-02c (#59 third PR, with N-06 #63):** `app.calibrate(acc, M | refArc)` over
  `timeSeriesPreparation.run` and `app.fitCalibration`; the 0.12.0 matrix interface is gone; A.2 reworded;
  referenceSystem >= 2.3. ☑ (merge run 38704980955: 414/414, predicted 414).
- **2026-10-15, N-02b (#59, second of three PRs):** `app.fitCalibration`, the estimator: CRN-filtered bias basis,
  model and K2 columns per axis through `lsqvce.fit` (kind A), VIF gate, block bootstrap, EIV bound; lsqvce in
  `dependencyReport`; 60 units. ☑ (414/414 on f9587d4, predicted 414, after two fixes the first run found).
- **2026-10-14, N-02a (#59, first of three PRs):** `app.syntheticPair` (the calibration truth through the real
  CRN chain) and `app.Calibration` (scale + bias curve, apply/invert/save/load); 59 units. ☑ (run 38685265825:
  408/408, predicted 408).
- **2026-10-14, N-01 (#58):** the six `qc` units forward to `timeSeriesPreparation.qc` with a
  `deprecated` warning; O-1/O-18 moved there; O-47 deletes them after the next release. ☑ with N-04 (run 38636885801: 394/394, predicted 394).
- **2026-10-14, N-04 (#61):** `sensor.Electrostatic` filters with the L1B CRN kernel from
  `timeSeriesPreparation.crnKernel("ACC")` by default, on 10 Hz dynamics whose grid it checks from
  `dyn.epoch`; `'OutputRate', 1` samples like ACT1B; `[]` stays unfiltered. timeSeriesPreparation and
  lsqvce are gate siblings since #76. ☑ (run 38636885801: 394/394, predicted 394).
- **2026-09-22, N-07 (#64):** `geom.grace` (GRACE-A/B, from R-5 = Bettadpur 2012 Sec. I.5.3) delivered
  on Tier B and ☑ on the CI runner (runs 35771831961 and 35773012179: 388 passed, 0 failed, 0 incomplete, predicted 388; the second's log is `tests/last_run.txt`), 57 units. O-46 opened (classic MAS1B is spacecraft
  mass; R-5 prints no inertia). The other eleven N-issues wait on timeSeriesPreparation, lsqvce#110 and N-02.

- **379/379 GREEN on Tier A** (2026-09-17, LP-E complete, R2026a Update 4): `matlab -batch ciGate`
  over an isolated dependency root, predicted 379, 0 failed, 0 incomplete, 0 category shortfalls,
  Code Analyzer clean on `+NonGravForces`, `auditPackage.py` clean (56 units, 56 catalogue rows,
  version 0.12.0). Rounds through the level: 374 (WP-E1, E3, E7), 379 (E2), 379 (E5, E6).
- **LP-E closed five obligations and opened two.** O-16 (SESAM) closed by implementing it; O-41
  (the GROOPS licensing question) closed the day it opened, because the data is public and is
  fetched rather than vendored; O-42 closed in two halves, the preset by owner ruling and the
  extrapolation warning by a second; O-38 and O-39 earlier. Opened: O-43 (the local rehearsal
  rewrites the developer's MATLAB path) and O-44 (`ceresClimatology` has no native evaluator and
  is refused at construction).
- **What LP-E added to the physics**: `io.importGroopsGrid` (both GROOPS containers, decoded from
  the bytes); `erp.backend = "ceresClimatology"`, the 1 degree CERES monthly climatology on the cap
  tiles; `drag.accommodation = "sesam"`, on both the MATLAB and the C path; and a warning when a
  computed accommodation leaves the range it was fitted in. `drag.gsi = "cll"` stays refused: O-17
  does not close, because [W14a]'s closed forms are whole-body and normal-incidence only and this
  kernel needs per-panel coefficients at arbitrary incidence.
- **THE PATTERN THIS LEVEL FOUND, and the reason to read Part G before writing anything**: seven
  defects were the same shape -- **a fact recorded in one place and compared with nothing**.
  `CLAUDE.md`'s unit count; the theory document's "51 units" and its `Model.rhs`; a data fetch
  printing `failed` three lines above `GATE TRUE`; the dependency minimum; the citation
  "JSWSC 9, A36", which resolves to no paper; `modelReference`'s "measured median 7.9e-9" beside a
  run file recording 5.475e-9; and `info.insideFit`, computed since 2026-09-07 and read by nobody.
  Each now has a mechanical check or a caller. **When looking for the next defect, look for a
  number that is written down twice and compared never.**

## 5. Superseded status (2026-09-16, after LP-C)

- **365/365 GREEN ON BOTH TIERS** (2026-09-16, LP-C complete, R2026a Update 4): on Tier A by
  `matlab -batch ciGate` over an isolated dependency root, and **on the CI runner** (run
  `35054822337`, both jobs, `GATE TRUE`). Predicted 365, 0 failed, 0 incomplete, 0 category shortfalls;
  Code Analyzer clean on `+NonGravForces`, `auditPackage.py` clean (55 units, 55 catalogue rows,
  version 0.10.0), `mex/test/runCChecks.sh` clean (14 bounds). Earlier rounds of LP-C: 334
  (C1–C4, C7), 352 (C8, C10), 358 (C11). **`Model.rhs` is gone** and the integrator's entry point
  is `app.DynamicsSource`, whose native C evaluator agrees with `Model.accel` to 2.5e-14 relative
  at a table node in six configurations.
- Before LP-C: **276/276 GREEN on the CI runner** (2026-09-08, LP-7 WP-30, R2026a Update 5, 219.6 s; predicted
  276). Before that LP-6b at 270/270, and: **262/262 tests GREEN** (2026-09-08, R2026a Update 4, Tier A host run; predicted 263,
  261 on the branch and 262 after merging main — reconciled in HANDOVER Part I),
  Code Analyzer 0 findings, gate TRUE;
  log in `tests/last_run.txt`.
  0.5.3 (the O-12 discharge) was written on Tier B and verified separately on the runner at
  245/245, exactly its own prediction; its test `tNgfIo/readACC1BMultiEmptyRobustness` is in
  this run too. The debt sweep was 244, LP-5x 238, LP-5 231, LP-4 169, LP-3 142,
  LP-2 116 (2026-09-06), LP-0 74 (2026-08-27).
  The six new tests came from the 2026-09-07 debt sweep: the gate's shortfall filter had
  been excusing any unit that held one exemption, hiding four missing category slots.
- Implemented (☑): L0 core, `geom.MacroModel`, `geom.graceFO`, `qc.hampel`,
  `qc.eventTable`, `sensor.Sensor` interface, test infrastructure, CI, and the
  LP-2 data I/O: `io.readACC1B`, `readSCA1B`, `readGNV1B`, `readTHR1B`,
  `readMAS1B`, `massModel`, `spaceWeather` (tests in `tests/tNgfIo.m`, real-data
  fixtures in `data/samples` and `data/spaceweather`), and the LP-3 environment:
  `env.density` (exponential, msis00 through the MEX of the vendored Brodowski C in
  `mex/`, built by `NonGravForces.buildMex`), `eclipseFactor` (cylindrical, conical,
  lunar), `earthRadiance` (pointSource, knocke cap tiles), `radiationSources`;
  pure-Python oracles in `tools/pyxcheck` (tests in `tests/tNgfEnv.m`); `eclipseFactor`
  solaarsCF (Robertson 2015, coefficients in `tools/pyxcheck/solaars_cf_coefficients.csv`);
  and the LP-4 geometry: `MacroModel` plate polygons and heat capacities, `graceFO` closed
  polygons, `boxWing`, `selfShadow` (ray cast + icosphere LUT), `AttitudeProvider` with
  `nominal`, `checkAttitudeSense`, `thermalState` (tests in `tests/tNgfGeometry.m`);
  and the LP-5 kernel and quality control: `kernel.gsiCoefficients` (sentman, maxwell,
  constantCd), `accommodation` (fixed), `dragPanelForce` (per species), `radiationPanelForce`,
  `thermalRecoil`, `assemble`, `bodyTorques`; `qc.residuals`, `detectImpulses`, `detectJumps`
  (PELT), `screen` (tests in `tests/tNgfKernel.m`, `tests/tNgfQcScreen.m`; oracles
  `tools/pyxcheck/sentman_check.py`, `pelt_check.py`); and the LP-5x cross-validation tooling:
  `tools/pyxcheck/ngfkernel/`, an independent standard-library Python mirror of the frozen
  kernel with its parity driver and `tests/tNgfMirror.m` (WP-33; WP-34 needs GROOPS, O-19).
- **No stubs left.** `io.importCERES` (WP-C8) and `env.magneticField` (WP-C10) were the last two;
  every catalogued unit is `done` except the `sensor.Sensor` interface row. The catalogue's
  `status` column is the authority (`doc/functionlist.csv`), not this sentence.
- LP-6 added five of the `Model` evaluation methods (`accel`, `torque`, `series`, `partials`,
  `budget`), `env.wind` (co-rotation) and `kernel.partials` (analytic derivatives). A sixth,
  `rhs`, existed from LP-6 to LP-C and is **removed**: the integrator's entry point is
  `app.DynamicsSource` (O-4, issue #35).
- Repository `github.com/astrobayer/NonGravForces` (the account was renamed
  from `matweigelt` on 2026-09-06; `referenceSystem` and `gracefoReading` moved
  with it), local clone at
  `C:\Users\matth\Documents\MATLAB\NonGravForces` with the dependencies as
  siblings.

---

## 6. Open tasks

**Done 2026-09-06: LP-2 (Data I/O)** except `readMAG1B` (no MAG1B file on the
host, HANDOVER O-8). Real-data fixtures: the five 2020-01-01 GRACE-C products
(1 Hz products are 600-record excerpts) and space-weather excerpts
2019-09..2020-03.

**Done 2026-09-07: LP-3 (env)** except solaarsCF (O-13), ceres (WP-27), msis21/
dtm2020/jb2008 (LP-3b), and the furnished-kernel Sun path (O-15).

**Done 2026-09-07: LP-4 (geom)** and solaarsCF (O-13 closed).

**Done 2026-09-07: LP-5 (kernel + qc)** except `partials` (LP-6), SESAM (O-16, **closed 2026-09-17 in LP-E WP-E3**),
CLL/Schamberg (O-17); O-2 closed (unit normals). Upstream issues filed:
referenceSystem#15 (O-10), gracefoReading#6 (O-12 — **fixed and closed upstream 2026-09-08**,
gracefoReading 1.0.1; the `l1bLoad` split is gone).

**Done 2026-09-07: LP-5x (WP-33)** — the kernel is mirrored by an independently derived
standard-library Python package and checked against it (`tests/tNgfMirror.m`). WP-34 needs a
GROOPS installation the host does not have (O-19).

**Done 2026-09-08: LP-6 (Model and partials)** — the whole chain evaluates end to end
against the real day, at 0.98 correlation with ACT1B cross-track. `rhs()` stays ◐ against
`gjToolbox`, which is not on the host (O-4).

**Done 2026-09-08: LP-6b (single-epoch C/MEX)** — ☑. Written on Tier B and verified through the
CI `matlab` job when the MCP bridge would not answer: predicted 270, ran 270, **270/270 GREEN** after
one 500 us absolute speed bound this session had set by guess was decomposed (not widened).
The 1e-14 parity is measured: 4.5e-16 on the force, 1.9e-15 on the torque, bitwise on `gamma`
and the speed ratios. The compiled kernel is 901.8 us per epoch against 4149.6 us for the MATLAB
chain, ratio 0.217 — for a 3-source scene, so it does not settle O-31. Both scope rulings were APPROVED by the owner on
2026-09-08: `thermalRecoil` is ported although the prompt did not name it, and `rsApiRotation`
is NOT linked because it is ITRS→GCRS while the kernel is body frame. D.6's speed justification was contradicted by O-29
and has been rewritten; under P-11 the speed ratio decides nothing and O-31 is informational.

**Opened 2026-09-08: LP-7** — `app.forceBudget` (WP-30) is delivered and GREEN (276/276 on the
CI runner, predicted 276): a chunked, multi-model
reduction over `Model.accel` whose oracle is `Model.budget` bitwise at one chunk. The two
dependency surfaces were read at last and both findings are in HANDOVER Part B: odeSuite is at
L0 with no dynamics interface and a signature that does not match `Model.rhs()` (O-4 narrowed),
and lsqvce 5.3's `fit` takes three declaration kinds with regularisation declared in the data,
which fixes `app.calibrate` as a kind-C batched declaration (O-24 discharged as a blocker).

`sensor.Electrostatic` (WP-29) is delivered too: bias, 3x3 scale, proof-mass offset coupling and
an FIR filter, with the omega sign convention pinned in its own help (`skew(omega) = R' Rdot` for
R body→inertial) because `referenceSystem`'s unit takes R the other way and carries a configurable
sense — the relation between them is RECORDED from the run, not guessed (O-32).

`app.densityInversion` (WP-32) is delivered too: the inversion is a DIVISION, because the gas force
is exactly linear in the species densities, so one model evaluation gives the whole sensitivity.
It does not calibrate and says so; it subtracts the modelled radiation and reports what is left
across the drag direction, which is the only honest quality flag it can offer.

`app.calibrate` (WP-19) is delivered with the TOOLBOX'S OWN SOLVER: the owner ruled lsqvce overkill
on 2026-09-09, and it is — the problem is linear and separates by axis, so it is three QRs of a few
dozen unknowns on an equilibrated design matrix. lsqvce is gone from the dependency table; O-24 and
O-33 are void.

`sensor.AtomInterferometer` (WP-29) is delivered: the triangular sensitivity kernel (a phase is an
average, not a sample), one reading per cycle with nothing interpolated across the dead time, and
the 2 pi ambiguity kept explicit — full phase, wrapped phase, and the fringe integer between them.
The Coriolis phase is driven by the ATOM velocity relative to the instrument (`dyn.atomVelocity`,
mm/s), not the orbital velocity, which atoms and satellite share: at 1.1 mrad/s and 1 mm/s it is
8.9 rad, an apparent 2.2e-6 m/s^2, seven times the drag. The 16.7 m/s^2 recorded here until
2026-09-15 came from feeding the orbital velocity in (LP-C F6).

`sensor.DragFree` (WP-29) closes the level: the residual a loop FAILS to remove, as a cascade of
bilinear high-pass sections with actuator saturation and a masked start-up transient.

**LP-7 IS COMPLETE; LP-R is open** — the C2 toolbox scan (`auditPackage.py` check 9) is delivered:
no call into Signal Processing, Statistics, Optimization, Curve Fitting, Symbolic, Aerospace or
Mapping, in code **or help text**, since an example that needs a forbidden product misleads the
reader who copies it. A name the file assigns to is treated as a variable, not a call.
**`tools/buildDocs.py` now generates every `doc/html` page and `helptoc.xml`**, each fact from its
own authority, and the audit refuses a stale one by CONTENT. O-11 is closed. It had already
failed: `Sensors.html` credited `Electrostatic` with a CRN filter the class deliberately does not
build, and called all three sensors future work months after two shipped.
**`demos/` is delivered and the gate RUNS every demo** (`tests/tNgfDemos.m`, discovered from the
folder): a demo nobody executes rots the same way that page did. Completion is asserted, not the
numbers — those are already bounded by the tests of the units the demo calls.
`tools/benchmark.m` regenerates the per-epoch costs (including O-29's 88 ms, which had been
measured once by hand and repeated ever since -- **and which was wrong: WP-C3 had already fixed it
to 1.21 ms on 2026-09-16 and nobody closed the row. Measured 2026-09-18, standard preset: 15.9 ms
per epoch at K=1, 0.798 at K=100, and a day at 1 Hz in 1 min 16 s via `Model.series`. O-29 and O-31
are closed; the table is committed at `tools/benchmark_reference.csv` and audit check 13 compares
every quoted figure against it**) and is deliberately NOT gated; `tools/packageToolbox.m`
builds the `.mltbx` from a programmatic options object and commits no `.prj`. Both live in `tools/`:
a development instrument is not part of the scientific surface. **O-35 answered** (2026-09-09): data is optionally downloaded during setup (`'FetchData'`, the same
three settings as `'FetchMsis'`) and otherwise the user supplies it. IERS, CelesTrak **and GRACE-FO
Level-1B** are all fetched — the last from the GFZ ISDC archive the owner supplied after my probe
guessed two wrong hosts and concluded it needed a login. **Two failed guesses at a URL are not a
finding.** The 2020-01-01 archive holds all 65 products, including MAG1B and ACT1B, which unblocks
O-8. `.mltbx` packaging is deprecated (O-36): publication is via git only. The four remaining stubs are blocked on data the project does not hold, so
they are not level debt. **Next is LP-R (release).** Two owner rulings of 2026-09-09 shape it: the
toolbox is **MIT** (`LICENSE`), and **O-29 was accepted** — 88 ms per epoch (superseded; see below)
in the standard preset shipped as a stated cost, not a blocker. That acceptance is moot: WP-C3 fixed
the cost on 2026-09-16 and the row went unclosed until 2026-09-18. **The standard preset costs
0.798 ms per epoch at K=100 and 15.9 ms at K=1**, the difference being a ~15 ms fixed cost per CALL;
`tools/benchmark_reference.csv` is the authority for that. **O-34 is closed** (2026-09-09): NRLMSISE-00 is no longer in this repository. `setupNonGravForces` clones it at
the pinned commit `bc9a2fe` when it is missing and owns the URL and the pin; `buildMex` resolves
`'MsisSource'` → `NGF_NRLMSISE00` → `mex/nrlmsise00` and REFUSES rather than fetching. **An
installer may reach the network, a build step may not.** The self-hosted Windows runner has no
bash, which is how the CI-step version died and how the better design was found. The oracle travels with the source (DOCUMENTATION §7), so tracking
upstream cannot silently invalidate a reference value — and that table turned out byte-identical
to the old copy's, so the move disturbed nothing.

**Done 2026-09-16: LP-C (completion), all twelve work packages, 365/365 on Tier A.** It answered a
whole-repository review's fifteen findings (F1–F15, HANDOVER Part G) and finished the `precise`
preset. What a later session most needs to know:
- **Self-shadowing was stated wrongly, not merely coarsely** (F1, F2). The table's "f = 0 on unlit
  plates" was a radiation-beam rule applied to a thermal gas; it deleted Sentman's edge-on shear on
  ~15 m² of plates and zeroed the Earth heating of every panel facing away from the Sun. New unit
  `geom.flowShadow` weights the occlusion by the incident drifting-Maxwellian flux. Measured on 120
  real epochs: drag moves 1.68e-9 m/s² RMS, thermal 1.77e-9, and the nadir panel stops reading 0 K.
- **`precise` now evaluates**: NRLMSIS 2.1, HWM14, CERES EBAF 4.2.1, `walker2014` accommodation (**since
  2026-09-17 it takes the FITTED `fixed` 0.85 of March et al. 2021 instead -- O-42**),
  IGRF-14. All third-party model source and data are FETCHED by `setupNonGravForces` at pinned
  SHA-256s and redistributed by nobody.
- **The integrator's entry exists** (WP-C9): `app.DynamicsSource` and `mex/ngfRhs.c`. `Model.rhs`
  is gone; `odeSuite` is a gate dependency. See §2's P11 and the dependency table.
- **Two instruments were found claiming checks they did not make**: the category audit credited a
  test to any unit whose name appeared in the suite name (F8; four collisions, so a unit could lose
  all its own tests and stay "complete"), and `CLAUDE.md`'s unit count was compared by nothing.
  Both now hold. **An instrument's silence is evidence only about what it inspects.**
- Open after it: O-38 (the theory PDF is behind its `.tex`; no LaTeX on this host), O-3 (CERES
  SYN1deg), O-17 (CLL/Schamberg — its closed forms are whole-body and normal-incidence only, so they
  do not serve a per-panel kernel), O-19 (GROOPS), O-26 (tagging). **O-16 closed 2026-09-17**:
  `drag.accommodation = "sesam"` ships on both the MATLAB and the C path.

**Done 2026-09-17: LP-E (sources), seven work packages, 379/379 on Tier A.** It took the gas-surface
and Earth-radiation models to their primary sources, which mostly meant correcting what this
repository had written down about them.
- **`io.importGroopsGrid`** reads both GROOPS containers, decoded from the bytes with no spec in
  hand. Its latitudes are **GEOCENTRIC**, and its help said geodetic until the `ceresClimatology`
  backend measured the file's own radii and found 0.168 m of spread against geocentric and 71.6 m
  against geodetic. A consumer that converts before the lookup misplaces every tile by 0.19 deg.
- **`erp.backend = "ceresClimatology"`**: the 1 degree CERES monthly climatology, measured exitance
  rather than an emissivity, piecewise constant in time. Refused at construction on the native path
  (O-44) rather than letting `find()` return empty and pack `none`.
- **`drag.accommodation = "sesam"`** closes O-16, on both paths, against an oracle that INTEGRATES
  eq. (26a) where the toolbox evaluates the closed form (26b). SESAM's pressure is a RAM pressure,
  140x the `n_O kB T` its neighbours use; feeding it theirs would have been silently wrong.
- **O-42 closed in two owner rulings**: all three presets take March et al. (2021)'s fitted
  alpha = 0.85, and a computed accommodation outside its fitted range now WARNS. It fires on
  GRACE-FO, which is the point -- the adsorption models give 0.11-0.17 at 520 km where the fitted
  literature reports 0.82-0.93, and they reproduce 0.86-0.95 at CHAMP's altitude, where the
  literature was fitted. The models are right; the extrapolation is what needed saying.
- **O-17 does NOT close**: [W14a]'s closed forms are whole-body and normal-incidence only.
- **O-6 narrows and stays open.** [H24] and [V] were both expected to supply a TRP oracle and
  neither does; what they gave instead is a magnitude check this toolbox passes (5.263 nm/s^2
  cross-track, between 4.22 for GRACE-FO and 6.4 for GRACE-A) and a published DISAGREEMENT about
  whether instantaneous reradiation helps at all (HANDOVER Part G). GROOPS `SimulateAccelerometer`
  is the standing ask.
- **The pattern**: seven defects of one shape, a fact written down twice and compared never. See
  section 5. **Next is LP-R (release).**

**LP-5x cross-validation** (owner rulings of 2026-08-27):
- WP-33: standard-library mirror of the **kernel layer only** under
  `tools/pyxcheck/ngfkernel/`, written from the equations (not transcribed from
  the MATLAB), **frozen once LP-5 is GREEN** as a regression oracle.
- WP-34: GROOPS as a **black box** (GPLv3 — no source transcribed in any
  language). `SimulateAccelerometer` output for the real day is frozen as CSV
  under `tests/fixtures/groops/` and committed; neither CI nor the host needs a
  GROOPS installation.

**Obligations**: `HANDOVER.md` Part H holds them, currently O-1 to O-50 with an owner
and a re-entry condition each. This file deliberately does NOT list them: the copy that
used to live here had drifted to O-15 of 21 and showed two closed items as open, and no
check compared the two (R1). Read Part H.

Housekeeping: **both halves of this paragraph are stale and kept as a warning.** The repository has
carried `LICENSE` (MIT) since 2026-09-09, and NRLMSISE-00 is no longer vendored at all (O-34
closed the same day) — `setupNonGravForces` clones it at a pinned commit. Every third-party model
this toolbox uses is now fetched, never redistributed, and each is recorded in
`THIRD-PARTY-NOTICES.md`. Corrected 2026-09-16 (LP-C WP-C12); it had been wrong for a week because
nothing compares this file with the repository.

**Version rule (owner's ruling, 2026-09-21 — supersedes the local rule of
2026-09-06):** the numeric fields of `toolboxVersion` are the authority and
`semantic` is composed from them. **MAJOR is a breaking interface, MINOR is
added functionality, PATCH is a fix** — one scheme for every repository,
defined in ci-actions `VERSIONING.md`, which is the authority this file
merely states. The local rule it replaced (minor = levels delivered, patch =
HANDOVER amendment) is dead; 0.12.0 is not renumbered, and the meaning
changes at the next bump. Note that this repository's one DELIBERATE, argued
deviation from `VERSIONING.md` is a different clause and is untouched: the
`gracefoReading >= 1.0.1` PATCH minimum of O-12/O-25.
`auditPackage.py` compares every copy
(Contents.m, README, doc/html footers) and checks `FunctionList.html` against
`tools/buildDocs.py`, which regenerates every page of `doc/html` and `helptoc.xml`; run both before every push.

---

## 7. Known pitfalls

Every one of these was an actual defect found in a predecessor or a dependency.

- **MAS1B is cold-gas tank mass** (~15 kg per tank, ~30.5 kg total), **not
  spacecraft mass**. Use the handbook constant: 601.214 kg (GRACE-C),
  601.209 kg (GRACE-D). Misuse gives 20–40× accelerations. `io.massModel`
  defaults to the constant and asserts > 500 kg.
  **Classic GRACE is the opposite**: its MAS1B (RL02) carries the SPACECRAFT mass,
  476.85 kg for GRACE-A on 2008-01-01 against the 487.2 kg launch mass of R-5
  Sec. I.5.3.3, so `io.massModel` cannot serve it (O-46).
- **R-5's surface table IS D-56935's**: the GRACE Product Specification (Bettadpur 2012)
  and the GRACE-FO handbook print the same eight rows digit for digit, and
  `tNgfGeom/graceMetamorphicAcrossDocuments` compares the two transcriptions. The
  dimensions differ (3123/720/694 mm vs 3122.5/721/695) and so does the origin: R-5's
  Fig 3 puts the C.o.G. 275 mm above the nadir plane, not at mid-height.
- **GNV1B is Earth-fixed (ITRS)**, not inertial. Rotate via
  `referenceSystem.rotationItrsToGcrs` + `applyRotation` (which also gives the
  velocity). A GPS tag is not UT1.
- **SCA1B is ICRF→SRF**; conjugate exactly **once** to get BODY→GCRS. Check
  against `referenceSystem.attitudeConvention`, do not assume.
- **Sentman**: name the coefficient before quoting a power of s. In this toolbox's
  normal/tangential convention the TANGENTIAL coefficient is l·[E/(s√π) + γZ] and has
  **no 1/s² term**; the lone Z/(2s²) belongs to the LIFT coefficient of the drag/lift
  pair. (v1 divided its shear by s instead of s² in the drag/lift convention it used,
  ×8 wrong at LEO s≈8 — HANDOVER C.1 D1.) The ambient gas is atomic O, not air:
  `R = 287 J/kg/K` is wrong by √(29/16) = 1.35 at 500 km.
- **One solar constant**: `const.S0 = 1361` W/m² everywhere. 1367 (Behzadpour
  Table 1) is only an `Options.srp.tsi` override for legacy comparisons.
- **Geodetic conversion must use the ellipsoid**; a spherical Earth gives
  −21.4 km polar altitude ⇒ density ×1.40.
- **CERES seasonal phase** references **Dec 22**, not Jan 1 (an off-by-10-day
  bug in an earlier toolbox).
- **GCRS right ascension ≠ sidereal time**; they differ by precession (~0.34°
  in 2026). Reference tests need the classical equinox route.
- **`exist()` does not resolve package-qualified function names.**
- **Quality, character and `_id` columns belong in `aux`, not `obs`**
  (`gracefoReading` convention).
- `gracefoReading.writeGroopsInstrument` / `readGroopsInstrument` were never
  verified against a real GROOPS installation — WP-34 closes this.
- `clear classes` wipes workspace variables inside a batched MATLAB script;
  split the script if state must survive.
- **Every oracle generator must run on both tiers**: no numpy anywhere (the host has none
  and the ubuntu audit job installs nothing). `hampel_check.py` had imported it and could
  not be regenerated at all until 2026-09-07.
- **`case` is not usable as a MATLAB table variable name**; the mirror fixtures use `caseId`.
- **The full suite exceeds the MATLAB MCP call budget**: run `matlab -batch ciGate` from a
  background shell and read the log; use the MCP session for single classes.
- **`referenceSystem` trims EOP rows before 1972-01-01 and warns** (since 2.1.216,
  2026-09-07; it refused outright before, O-10). The whole `eopc04_20.1962-now` may
  now be passed to `readGNV1B`, with 3652 of 23562 rows ignored.
- **MAS1B samples start at 3599.99 s of the day**, so the first hour has no
  sample; `massModel` holds the nearest value for two sample intervals.
- **Fixture excerpts are real data with `num_records` rewritten**; a test that
  skips for a missing file counts as `incomplete` and reddens the gate, so
  fixtures are committed, never fetched.
- **A fixed global lat-lon grid for Earth radiation misses the flux by half at
  490 km** (the sub-satellite tile dominates 1/rho^2); `earthRadiance` tiles the
  visible cap in rings around the sub-satellite point (`erp.grid` = [rings sectors]).
- **The Montenbruck-Gill lens area with `acos` loses 1e-11 in nu** when the Earth's
  disc is 250x the Sun's; `eclipseFactor` uses atan2 arcs (tested to 1e-12).
- **The published NRLMSISE-00 test cases use lst = 16 h regardless of UT and
  longitude**; reproduce them only with `density(..., 'LocalSolarTime_h', 16)`.
- **MATLAB forgets `MW_MINGW64_LOC`**; `NonGravForces.buildMex` locates the MinGW
  support package itself and runs `mex -setup C` for the session.
- **Angles by atan2(|a×b|, a·b), never acos of a dot product**: acos loses sqrt(eps)
  near 0 and 180 degrees (1e-6 deg); bitten twice (eclipse lens area, attitude sense).
- **The handbook's exterior dimensions do not close at the 2-4 cm level**; `graceFO`
  derives the nadir width from the 50 degree flank crossing so the polygons do not
  intersect (a spurious 2.5 % self-shadow from straight above was the symptom).
- **A BLAS matrix product sums in a batch-size-dependent order**; use column-wise sums
  where a vectorisation test must be bitwise.
- **The owner's host has Signal Processing and Statistics; the CI runner does not.** Anything
  from them (`corr`, `filtfilt`, `prctile`, ...) passes locally and fails on CI. `auditPackage.py`
  scans for them since 2026-09-08; base MATLAB and the PCT only (HANDOVER D.7 C2).
- **`max`, `min` and `sum` omit NaN by default**, so a test bound can read 0 while the unit
  returned NaN. `NgfTestCase.verifyBoundAndRecord` now refuses a non-finite measurement; that
  is how a four-level-old NaN in `geom.selfShadow` was finally seen (2026-09-08).
- **Sentman's shear is lift**: the gas force is NOT antiparallel to the relative velocity
  (0.34 deg off on GRACE-FO). Only `constantCd` is, by construction.
- **The co-rotation term is `v - Rdot*(R'*r)`**, from `env.wind`, never a hand-written omega x r.
- **Sentman's published C_D/C_L are drag and lift, not normal and tangential**; the kernel
  derives Cp (normal) and Ctau (tangential) from the momentum integrals and asserts both
  projections. Three independent derivations reached the same conclusion in LP-5x. The limit
  that settles it: edge-on, Cp = 1/(2s^2) and Ctau = l/(s*sqrt(pi)). The theory document §2.3
  mislabelled them until 0.10.0; `doc/src/*.tex` is corrected and **the committed PDF is not
  rebuilt** — no LaTeX on the host (O-38).
- **The speed ratio is per SPECIES, not from a mean molar mass**: the mean-mass form is
  harmless on a ram face, -24 % at 100 deg incidence and collapses at 135 deg (measured 2026-09-07).
- **1 + erf(x) cancels for x << 0**: use erfc(-x). **1 + a(x - 1) loses 1e-8 at a = 1**: write
  (1 - a) + a x.
- **Hampel flags the peaks of a strong pulse** (runs of several samples); replace only
  isolated runs (<= 2 samples) before a matched filter, and let pulses claim their samples.
- **BIC for change points counts the location**: (p + 1) sigma^2 log K, not p sigma^2 log K;
  PELT must keep not-yet-admissible candidates in its set.
- **`MacroModel.N` holds unit normals** (renormalised); `panels(k).normal` keeps the printed
  handbook digits (norm 1 - 7.3e-7, measured).
- **`M_PI` is not in C99** and is absent under MSVC without `_USE_MATH_DEFINES`; `ngfKernel.c`
  defines its own `NGF_PI`.
- **A SCALAR indexed by an empty index is 1x0, a column indexed by it is 0x1.** `x(sIdx)` with `sIdx = P + (1:0)`
  breaks `J * x(sIdx)` only when x happens to be a scalar -- in `fitCalibration`, a constant bias on an axis
  with no estimated scale. Shipped in N-03a, found by N-03b. Force the shape: `reshape(x(sIdx), [], 1)`.
- **A dead housekeeping channel need not be constant**: AHK1B's TICUR reads -273.2 C with a 3 mK wobble, so a
  constancy check passes it. Look at the values before choosing a regressor (N-11a).
- **lsqvce's `B` spans every local column of the system**, not just the arc's own day: m_a x nl over all days,
  the day's block in place (sparse). And a system needs at least one COMMON column (N-09).
- **Vectorising is not by itself faster.** `radiationPanelForce` broadcast to one 3xNxK expression per source
  was 18 % SLOWER than the plate loop it replaced; one NxK accumulator per component was 28 % faster. And the
  profiler showed the faster version as slower: compare variants interleaved with tic/toc, never by profile.
- **A table that does not converge with its grid is mis-indexed, not coarse.** The flow-shadow prototype's
  error stayed at 0.1 from 21 to 41 nodes because nodes and queries used mirrored coordinates. Refine once
  before believing any interpolation error.
- **A title's text object hides its handle from `findobj`**; a test reading back a Tagged title needs `findall`.
- **`tsFitAR`'s default exact method is too slow for anything drawn**: over a minute for 4000 samples with a
  gap. Fit Yule-Walker (`Method = "yw"`) on a gap-free run and filter (N-05).
- **`struct('a', {})` builds a 0x0 struct ARRAY**, not a scalar with an empty field; only cell
  values distribute. Assignment into an element of a struct array also needs the same fields in
  the same ORDER, which is why `buildMex/localState` rebuilds the record instead of extending it.
- **A MATLAB function call result cannot be indexed** (`f('dx')(1)` is a syntax error); a dynamic
  field can (`r.(name)(1)`), which is how the mirror fixtures are read.
- **`rsApiRotation` is ITRS-to-GCRS** and has nothing to do with the body-frame kernel; the plan
  said the single-epoch MEX would link it and the plan was wrong (HANDOVER B.1).
- **A dependency's surface is a guess until its source is read.** Two records in Part B described
  toolboxes that had moved on: lsqvce's flat 1.0 surface was gone, and gjToolbox has no dynamics
  interface at all. Both repositories are readable through the GitHub API from any tier.
- **The `fast` preset sets `thermal.model = "none"`**, so a thermal invariance asserted on a fast
  model asserts 0 == 0 and measures nothing.
- **`env.density("exponential")` returns species and T as NaN on purpose** — it models neither.
  Any unit that screens an input it does not read will refuse the fast preset's own environment;
  `kernel.singleEpoch` did exactly that until 2026-09-08. Build a test's inputs from the real
  units, not from literals, or this class of defect stays invisible.
- **Separate the algebra from the statistics, and put the tight bound on the algebra.** A claim
  about a covariance formula is an identity (`sigma_x * sqrt(K) / sigma0 = 1` for a constant
  column) and should be asserted at 1e-12 on one arc, not as a ratio between two noise
  realisations. Three runs were lost learning this.
- **1/sqrt(K) is a law for a CONSTANT design column and for nothing else.** A bias sigma falls
  that way; a scale sigma is set by the spread of the modelled signal over the arc, and a longer
  arc need not carry more spread. Asserting the law on the wrong parameter measured 8.3 against
  an expected 2.
- **A refusal test that does not refuse is not a passing test, it is an untested guard.** A 2 s
  knot interval was estimable and merely useless; only 0.5 s put the unknown count above the
  sample count. Provoke the branch you mean, and pair it with the call that must succeed.
- **Dropping a dependency leaves dead code behind it.** Removing the lsqvce row left its only
  helper unreferenced and `codeIssues` (DEFNU) caught it before any test ran. Grep for the
  helpers a removed row used, not only for the row.
- **An empty L1B day is a refusal, not an empty array**: gracefoReading raises
  `gracefoReading:readL1b:noData` for `num_records: 0` (GRACE-D 2020-01-19..21). A long run must catch that
  identifier and write each finished part as it goes -- the first month run lost three hours to it.
- **MATLAB's `websave` times out after 5 s by default.** A research server's first byte can take longer;
  `setupNonGravForces` passes `weboptions('Timeout', 60)` and retries once (O-50). Any new fetch goes through
  its `localDownload`, never a bare `websave`.
- **A status code is not a file.** `https://www.ngdc.noaa.gov/IAGA/vmod/igrf14.f` answers `200`
  with an HTML error page. Check the content, not the response code — met twice in two days.
- **Never claim bitwise across a change in BATCH SIZE.** Vectorised products and norms associate
  by batch size, so K = 1 and K = 10 are not the same arithmetic even when they are the same
  formula. Met three times now: `pagemtimes` in `Model.evaluate`, the chunked RMS in
  `forceBudget`, and one epoch vs an arc in `densityInversion`. Bitwise IS claimable across a
  change of grouping that leaves the arithmetic alone — two models stacked, a maximum of maxima.
- **A numerical bound written from an estimate is a red gate waiting to happen.** Three sessions
  in a row have lost a run to one. Where cancellation governs, COMPUTE the bound in the test from
  `eps` of the terms involved instead of writing a literal. Where the model is exactly linear, a
  LARGER finite-difference step is strictly better, the opposite of the usual rule.
- **The CRN kernel is defined on the 10 Hz grid.** 1407 taps at 0.1 s; convolving a 1 Hz arc with them is a
  filter of ten times the span that nothing would report. `sensor.Electrostatic` therefore requires
  `dyn.epoch` with its default kernel and refuses any other spacing; for unfiltered 1 Hz work pass `[]`.
  Its DC gain is 1 + 1.115e-7 (independent Python evaluation of ATBD App. B: 1 + 1.11516e-7), and it
  is normalised at f0, NOT at DC -- never renormalise it.
- **`skew(omega) = R' Rdot` is the sensor's convention**, and
  `referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference")` returns its
  NEGATIVE (ratio -1.0, measured). Reading the skew entries off `R' Rdot` without projecting the
  symmetric part out leaves the differencing error in (6e-7 relative).
- **A rehearsal that runs FEWER checks than the runner is not a rehearsal either.** `ciGate` does
  not run `codeIssues`; the CI `matlab` job runs it as a separate step BEFORE `ciGate`, so a local
  gate that only calls `ciGate` can be 379/379 green and still redden CI. #46 did, on one INUSD
  warning, after an argument stopped being used when a geodetic conversion was removed. Run
  `codeIssues(fullfile(pwd, "+NonGravForces"))` before believing a local gate. This is the mirror
  of the entry below, which is about giving the local run MORE than the runner gets.
- **A rehearsal that gives the local run a hint CI does not get is not a rehearsal.** `buildMex`'s
  walk to odeSuite's `csrc` overshot by one directory (`fileparts` x3 already reaches the repository
  root; the extra `'..'` left it). Four green gates in a row said nothing, because every local script
  set `NGF_GJTOOLBOX` and CI sets nothing — the branch CI depends on had never once run. Clear the
  environment variables the runner lacks before believing a local gate, and TEST the fallback.
  **The same trap has a second door: a build ARTEFACT the host happens to hold.** `gjToolbox`'s
  native stepper is not committed, so a fresh checkout has none, and `tNgfIntegration` propagates
  through it. Two green gates hid that too, because this host's odeSuite clone had been built by
  hand. `ciGate` now runs the dependency's own `buildtool mex` (never a recipe copied here), and the
  rehearsal deletes `gj_stepper` first. **Ask what the runner does NOT have that you do.**
- **A pin taken from the tool whose work it checks is not a check.** The DE440 digest was first
  copied out of `referenceSystem.downloadKernels`' own manifest; it matches no file on disk (O-40),
  so setup shipped a pin that failed on the very first CI run. Measure the artefact yourself.
- **Existence is not content, and a check that runs only on download protects only the first run.**
  `localFetchKernel` verified the SHA-256 when it fetched and returned "present" otherwise, so on
  every machine that already had the file the pin did nothing - proved by injecting a deliberately
  wrong pin and watching the gate stay green. It hashes on every call now.
- **A failed check must destroy what it rejected, and the gate must read the report.** On
  2026-09-16 CI printed `DATA de440 failed: SHA-256 ...` three lines above `GATE TRUE`: the
  unverified kernel stayed exactly where every consumer looks, the tests read it and passed, and
  nothing in `ciGate` had ever looked at `info.data.items`. Both are fixed and fault-injected.
- **referenceSystem's analytic Sun is good enough for radiation pressure, measured**: against JPL
  DE440 it is 0.0031 arcsec in direction and 0.3 km in range, i.e. **4.0e-9** of the solar radiation
  pressure - three orders below the 1e-3 uncertainty of the solar constant. `rs.de440` is the
  `precise` default and is now fetched and tested (O-7), but nobody should expect it to move a
  number. **It had never been evaluated by any test until 2026-09-16**, only its refusal had.
- **Never assert a host-dependent absolute** (rule T5). A 500 us per-epoch bound written from a
  guess reddened the gate on a runner that measured 901.8 us, in the same session that quoted the
  rule. Assert the ratio; record the absolute with a ceiling an order of magnitude above it.
- **When the MCP bridge is down, CI is still a MATLAB tier**: push the branch, let the `matlab`
  job run `ciGate`, read `tests/last_run.txt` from the run artefact. That is how LP-6b was verified.
- **ASK WHETHER THE PR MERGED BEFORE PUSHING A FOLLOW-UP TO ITS BRANCH.** `main` is not protected
  against this and git will not warn: when a PR is merged the branch is deleted, and a later `git
  push` to that name SILENTLY RECREATES it. The push succeeds, `git ls-remote` confirms the SHA, and
  **no CI runs at all**, because no open PR points at it — a green-looking push that was never
  tested. The tell is `* [new branch]` in the push output for a branch that already existed; read
  it. This happened twice on 2026-09-18, once after #54 merged and once after #55. Recovery is
  cheap (branch from HEAD, push, open a new PR, delete the recreated one) but only if it is noticed,
  and the failure mode is indistinguishable from success in every check except the run list.
- **RE-RUN `auditPackage.py` AFTER THE GATE, NOT BEFORE IT.** The gate writes `tests/last_run.txt`
  and check 11 compares `README.md`'s quoted test count against it, so an audit run before the gate
  is an audit of the previous run's numbers. On 2026-09-18 that reddened CI on a branch whose local
  audit had been clean twice: two tests were added, the gate recorded 381, `README.md` still said
  379, and nothing had re-read the pair. **A check is only evidence about the moment it ran.** The
  same mistake in a different costume cost a CI wait the same day: a loop that asked "are any checks
  pending?" was satisfied by the PREVIOUS commit's completed checks and reported green for code CI
  had never seen. Pin the query to the SHA, and re-run the cheap check last.
- **ONE WARM-UP CALL IS NOT WARM.** `tools/benchmark.m` warmed once before timing, and on 2026-09-18
  that cost a factor of 5.4 on `kernel.singleEpoch` (2.44 ms against a settled 0.456) and 7.8 on the
  fast preset at K=1 (31.9 against 4.00). The evidence had been visible for an hour without being
  read: three consecutive invocations in one MATLAB process, the first reporting the high value and
  the next two agreeing on the low one. It warms five times now. **An argument had already been
  built on the artefact and had to be withdrawn** — that the presets differ by 4 % at K=1, when they
  differ by 4.0x — which is the real cost of a bad measurement: not the wrong digit, the conclusion
  drawn from it and written somewhere else.
- **A per-epoch cost is mostly a property of the BATCH, not of the toolbox**, so never quote one
  without K. `Model.accel` in the standard preset costs 11.9 ms at K=1 and 0.440 ms at K=100
  (0.16.0, 2026-10-18, with the flow-shadow table), because one call pays ~11.6 ms of setup before any
  per-epoch work (fast: 5.6 ms and 0.063). A day at 1 Hz through `Model.series` was 1 min 16 s on 2026-09-18 -- a figure
  `tools/benchmark.m` does NOT produce, so nothing regenerates it. The authority is `tools/benchmark_reference.csv`,
  which `tools/benchmark.m` writes and `auditPackage.py` check 13 compares every quoted figure
  against; `tools/benchmark_last.txt` is the gitignored host log and is NOT the authority, because a
  check whose authority is absent in a fresh clone protects nothing.
- **An OPEN obligation row is a standing instruction to repeat what it says.** O-29 was discharged by
  WP-C3 on 2026-09-16 and left open; over the next two days its stale 88 ms was copied into four more
  files and one shipped PDF chapter, while the correct value sat in an executed test. **Closing a row
  is part of discharging it**, and a work package that fixes a number must be the thing that deletes
  the old one. Checking that a fact appears once does not catch this: the second copy was an
  obligation to change the first.
