# HANDOVER — NonGravForces 0.3.0 (concept 3)

**Single source of truth** for the design, the audit that motivated it, the work plan, the rules, and the decisions of `NonGravForces`, a MATLAB toolbox computing non-gravitational forces, torques and partial derivatives for LEO satellites, with GRACE-FO accelerometer processing as the reference application and a sensor layer that admits future instruments (electrostatic, cold-atom, optomechanical, drag-free).

Revision 2.2 (LP-3 Environment delivered; version 0.3.0) · 2026-09-07 · Claude Fable 5.1 (Anthropic) for Matthias Weigelt · MATLAB R2026a Update 4. Revision 2.1 (LP-2, 0.2.0) · 2026-09-06. Revision 2.0 · 2026-08-27 · Claude Fable 5 (Anthropic).

Companion: `doc/NonGravForces_Theory_and_Implementation.pdf` (theory, practical considerations, implementation). This document does not repeat the physics; it cites the PDF by section.

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
| Quality control: `qc.hampel` (with numpy oracle), `qc.eventTable` (schema) | `+NonGravForces/+qc/`, `tools/pyxcheck/` | ☑ |
| Sensor interface `sensor.Sensor` (abstract) | `+NonGravForces/+sensor/` | ☑ (interface) |
| Data I/O (LP-2, 2026-09-06): `io.readACC1B`, `readSCA1B`, `readGNV1B`, `readTHR1B`, `readMAS1B`, `massModel`, `spaceWeather`; real-data fixtures under `data/samples`, `data/spaceweather` | `+NonGravForces/+io/`, `data/`, `tests/tNgfIo.m` | ☑ (116/116, 2026-09-06) |
| Environment (LP-3, 2026-09-07): `env.density` (exponential, msis00 via the MEX of the vendored Brodowski C), `eclipseFactor` (cylindrical, conical, lunar), `earthRadiance` (pointSource, knocke cap tiles), `radiationSources`; `buildMex`; pure-Python oracles in `tools/pyxcheck` | `+NonGravForces/+env/`, `+NonGravForces/buildMex.m`, `mex/`, `tests/tNgfEnv.m` | ☑ (142/142, 2026-09-07) |
| 28 documented stubs raising `NonGravForces:<name>:notImplemented` | all packages | ○ |
| Test infrastructure: `NgfTestCase`, `runNgfTests`, `auditTestCategories`, exemptions; 7 test classes, 142 tests, **gate GREEN 142/142** (74/74 at LP-0, 116/116 at LP-2) | `tests/` | ☑ |
| CI: `ci.yml`, `ciGate.m`, `auditPackage.py`; Code Analyzer clean (0 findings) | repository root, `.github/workflows/` | ☑ (local run 2026-08-27) |
| Help-browser documentation (`info.xml`, `helptoc.xml`, six pages) | `doc/html/` | ☑ (`doc NonGravForces` opens) |
| Theory & implementation PDF | `doc/*.pdf` (+ LaTeX source in `doc/src`) | — |

Run log of the verifying run: `tests/last_run.txt` (2026-09-07, R2026a Update 4, 142 passed, 0 failed, 0 incomplete, 0 audit shortfalls, predicted 142; LP-2: 116/116; LP-0: 74/74).

### A.2 The rule that governs every other rule

> **NonGravForces computes forces and torques and nothing else.**

Four toolboxes of the same owner form the ecosystem; each is called by qualified name and none is copied, vendored, wrapped or aliased:

| Toolbox | Owns | NonGravForces uses |
|---|---|---|
| `referenceSystem` ≥ 2.1 (2.1.214 in hand) | time scales, EOP, ITRS↔GCRS with derivatives, geodetic coordinates, quaternion algebra, attitude interpolation/integration, Sun/Moon/planet states (SOFA analytic or DE440 via MICE) | everything kinematic |
| `gracefoReading` ≥ 2.1 (2.3.0 in hand) | L1B ASCII reading, quality flags, GROOPS export | all L1B input |
| `gjToolbox` (odeSuite) | Gauss-Jackson orbit integration; `gjToolbox.Dynamics` interface; force-evaluation budget | `Model.rhs()` returns the non-gravitational acceleration handle a caller composes with its gravity model and hands to the integrator |
| `lsqvce` 1.0 | batch least squares, VCE, ARMA whitening, data snooping, robust reweighting, diagnostics | `app.calibrate` builds design matrices and delegates the adjustment; `qc` supplies masks and jump columns |

Consequences that differ from v0.1.1: **no gravity, no equations of motion, no integrator, no least-squares code** inside NonGravForces (v0.1.1 planned WP-17 `dyn/gravity` + `eomFcn` and hand-rolled calibration; both are removed, see C.3).

### A.3 How to run

```matlab
run('C:\Users\matth\Documents\MATLAB\referenceSystem\setupReferenceSystem.m')
run('C:\Users\matth\Documents\MATLAB\gracefoReading\setupGraceFOReading.m')
cd('C:\Users\matth\Documents\MATLAB\NonGravForces'); setupNonGravForces   % adds root, tests/, tools/; registers doc
runNgfTests                          % gate; 'ExcludeSlow',false adds the speed tag
doc NonGravForces                    % help browser
```
`setupNonGravForces` adds `tests/` itself (v0.1.1 did not: finding V-13 friction, closed by design).

**Continuous integration** (`.github/workflows/ci.yml`, modelled on referenceSystem's): job `audit` runs `auditPackage.py` (stdlib only: catalogue ↔ files, help template, `fail` taxonomy, banned constructs); job `matlab` checks out `referenceSystem` and `gracefoReading` as siblings, runs `codeIssues` on the package (must be empty) and then `ciGate`, the one entry point that installs, runs the suite and judges. `ciGate` is also the local postflight (`matlab -batch ciGate`). Branch `main` is protected: both checks required, no force pushes, changes arrive by pull request.

### A.4 Session protocol (unchanged in spirit from v0.1.1, tightened)

1. One level = one session = one deliverable zip. Paste the standing instruction (E.1), the level prompt (Part G) in full, and this document.
2. Nothing is ☑ without a MATLAB run report read by the author. `runNgfTests` writes `tests/last_run.txt`; the tick rests on that file, not on the pass count.
3. Predict the test count before the run; reconcile after (P-2).
4. A test failure is decomposed before any bound moves (R6). In this very session two failures were test-authoring errors (a typed oracle sum, a 0-based index) and one was a bound set by guess (interpreted vs compiled speed ratio); each is recorded in the test's comment with the measured number.

---

## Part B — Dependency contracts (measured 2026-08-27)

### B.1 referenceSystem 2.1.214 — functions NonGravForces calls

| Need | Call | Convention |
|---|---|---|
| epoch from L1B | `daySecondToEpoch(Nx2,"GPS")` | epoch struct `.day .sec .scale .source` |
| scale conversion | `convertTimeScale(ep,"TT"|"UT1"|"UTC")` | UT1 needs EOP |
| ITRS→GCRS | `rotationItrsToGcrs(ep,'Eop',eop)` → `[R,Rdot,Rddot]` | 3×3×N pages, ITRS→GCRS, tidal corrections default on, samples EOP beyond the arc (±3 d margin). **The series must hold no row before 1972-01-01** (measured 2026-09-06: the shipped `eopc04_20.1962-now` read whole raises `referenceSystem:leapSecondTable:beforeUtcEpoch` from `interpolateEop`; `finals2000A.all`, or a C04 subset from MJD 41317 on, works; the two agree to 9 mm / 5e-6 m/s on the 2020-01-01 arc). Cost measured: 600 epochs with `Rdot` 0.14 s, without 0.024 s (O-10) |
| apply | `applyRotation(R,Rdot,r,v)` | carries ω×r; never formed by hand |
| geodetic | `cartesianToGeodetic(Nx3,'Ellipsoid',"WGS84")` | **rows in, default GRS80** → WGS84 passed explicitly |
| Sun/Moon | `geocentricBodyState("sun",ep)` / `Source="mice"` | ICRS-oriented 3×N; analytic tier RMS 8.3e3 m; DE440 tier via furnished kernels (de440s.bsp present) |
| attitude | `interpolateQuaternions(4×M,nodeEp,queryEp,'Method',"slerp"|"squad")`, `rotationFromQuaternion`, `multiplyQuaternions`, `invertQuaternion`, `attitudeConvention()` | scalar-first Hamilton, passive; SCA1B is ICRF→SRF, NonGravForces stores BODY→GCRS, conjugated once in `geom.AttitudeProvider` |
| angular rate | `angularVelocityFromRotationRate(R,Rdot)` | CoM-coupling term in `sensor.Electrostatic` |
| orbit frame | `rotationGcrsToRsw(r,v)` | budgets and plots |
| constants | `iersConstants()` (`.speedOfLight .earthRotationRate .gmEarth …`), `ellipsoidParameters` | no `au`: NonGravForces owns it in `const` |
| version | `toolboxVersion().version` | **struct, not string** — the first probe of this session was wrong and was corrected by measurement |
| C API | `mex/rsApi.h` `rsApiRotation` (0.56 µs) | for the single-epoch MEX (LP-6b) |

### B.2 gracefoReading 1.0.0 (measured 2026-09-06; was 2.3.0 on 2026-08-27, renumbered by its owner)

`readL1bMulti(files)` → `[time, obs, aux]` with `aux.products`, `aux.fileNames`, `aux.recordsPerFile` (no `product` field); `checkQuaternionNorm(obs, aux)` → `[ok Nx1, stats]`, default tolerance 1e-6, measured 9e-11 on the real SCA1B sample; `readL1b(file)` → `[time Nx2 GPS, obs, aux]`; `readL1bMulti`; `decodeQualFlags`; `checkQuaternionNorm`; `findGaps`; `version()`. Rules learned in v0.1.1 and kept: `obs` holds measurements only (flags, char columns, `_id` columns, time tags go to `aux`); `prod_flag` products have per-record columns (bit 0 rightmost); **MAS1B is cold-gas tank mass (~30 kg), not spacecraft mass (601.214 / 601.209 kg)** — misuse gives 20–40× accelerations; `io.massModel` defaults to the handbook constant and asserts > 500 kg.

### B.3 gjToolbox (odeSuite)

Rule R9 of that toolbox: no force modelling there; force models implement `gjToolbox.Dynamics`. `Model.rhs()` returns `@(t,r,v) a_ng` (3×1, GCRS, m/s²) plus a `budget` counter hook; composing with gravity is the caller's dynamics object. Not on the path on 2026-08-27 (optional dependency; `dependencyReport` reports it absent without error).

### B.4 lsqvce 1.0

`adjust` / `Accumulator` (flat), `armaWhiten`, `snoop` (iterative snooping with Cholesky downdate), `robust` (Huber/bisquare/IGG3), `varComp`, `diagnose`, `report`. `app.calibrate` supplies: observation vector (sensor observable), design blocks (bias spline, 3×3 scale, jump columns from `qc.eventTable`), a-priori ARMA noise from the sensor `noiseModel`, and the `qc` mask; it receives `Results`.

---

## Part C — The audit (both predecessors, in every little detail that matters)

### C.1 SatelliteNonGravToolbox v1 — verified defects (numbers re-derived, see PDF §2)

| # | Defect | Evidence in source | Magnitude |
|---|---|---|---|
| D1 | Sentman: exponential group of C_p and the whole shear term divided by s instead of s² | `satDragAccel.m` lines `Cp = ... / s_ratio`, `Ctau = ... / (s_ratio*sqrt(pi))` | shear ×8 at s≈8 (LEO) |
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

### D.2 Layers (a layer calls only lower layers plus the four dependencies)

```
L8 app      calibrate (lsqvce bridge), forceBudget, densityInversion
L7 qc       residuals, hampel, detectImpulses, detectJumps, eventTable, screen
L6 sensor   Sensor (abstract), Electrostatic, AtomInterferometer, DragFree
L5 Model    accel | torque | series | partials | rhs | budget   <- single point of access
L4 kernel   gsiCoefficients, accommodation, dragPanelForce, radiationPanelForce,
            thermalRecoil, assemble, bodyTorques, partials
L3 geom     MacroModel, graceFO, boxWing, selfShadow, thermalState, AttitudeProvider, checkAttitudeSense
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
app.calibrate ──▶ lsqvce.adjust(A = [bias spline | 3x3 scale | jump columns], whitening from sensor.noiseModel, snoop) ──▶ Results + provenance
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

Vectorise over K epochs; environment in bulk (density backend accepts K points; MEX for the MSIS inner loop is the sanctioned C use); `pagemtimes` for attitude; the panel loop runs over N≈8 with K-vectorised arithmetic; self-shadow LUT replaces per-epoch ray casting; `series` splits K into chunks of `Options.series.chunk` and runs them serially or under `parfor` (PCT feature-detected) with results required bit-identical either way (test category `vectorization` + a `parallel` tag). The single-epoch C/MEX (LP-6b) ports only the per-epoch kernel (GSI + radiation panel loop + assemble) and links `rsApiRotation`; it exists for `Model.rhs()` inside an integrator, where MATLAB call overhead (5–15 µs) dominates. `gpuArray` is not used (N=8 panels give nothing to a GPU; ray tracing is the only candidate and is pre-computed into the LUT).

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
| env | ☑ (frozen 2026-09-07) `[rho, species 8×K, T 2×K, info] = density(model, ep, lla 3×K, sw, 'Drivers', 'AnomalousOxygen', 'Exponential', 'LocalSolarTime_h')` (exponential, msis00 via MEX; msis21/dtm2020/jb2008 refuse `notImplemented`); `[nu, info] = eclipseFactor(r, r_sun, 'Model', "cylindrical"|"conical"|"solaarsCF", 'Moon')` (solaarsCF refuses, O-13); `tiles = earthRadiance(backend, ep, grid, r)` (pointSource, knocke cap rings; ceres refuses); `[src, info] = radiationSources(ep, r, q_body2gcrs, options, 'Ephemeris', 'Sun', 'Moon')` with `src.dir` 3×M×K body (satellite → source), `.P` M×K Pa, `.band`, `.nu`. ○ `wind(model, ep, r, v, lla, sw)`; `magneticField(ep, r_itrs) → B` | ☑ 4 / ○ 2 | 05/25, 07, 08/27, 23 |
| geom | `MacroModel` ☑, `graceFO` ☑, `boxWing(dims, optics)`, `selfShadow(macro, dir_body, 'Method', "raycast"|"lut")`, `thermalState(macro, src, 'Model', …) → T_panel`, `AttitudeProvider(ep_nodes, q_nodes)` with `.at(ep)` (slerp/squad via RS), `checkAttitudeSense(q, r, v)` | ○ | 09–11, 28 |
| kernel | `gsiCoefficients(model, s, theta, alpha, Tw, Tinf) → [Cp, Ctau]` (Sentman exact erf form, CLL, Maxwell, Schamberg); `accommodation(model, species, Tinf, …) → alpha`; `dragPanelForce(macro, rho, species, v_rel_body, …) → F_k 3×N×K`; `radiationPanelForce(macro, src) → F_k`; `thermalRecoil(macro, T_panel) → F_k`; `assemble(F_k, macro, r_cm) → [F, tau]`; `bodyTorques(macro, r, B, …)`; `partials(…)` analytic + complex-step verifier | ○ | 12, 13, 26, 24, 16 |
| sensor | `Sensor` ☑ (abstract: `observe`, `partials`, `noiseModel`, `numParameters`); `Electrostatic(bias, S, dr_proofMass, filter)`; `AtomInterferometer(k_eff, T, cycle, axis, 'Rotation', …)`; `DragFree(controller)` | ○ | 29 |
| qc | `hampel` ☑, `eventTable` ☑, `residuals(y_meas, y_model, sensor, cal0) → r_white`; `detectImpulses(r, 'Response', crn, 'Priors', thr) → events`; `detectJumps(r, 'Model', "level"|"slope", 'Priors', eclipse/mass) → events` (PELT + BIC); `screen(…) → [mask, events, info]` | ○ | 31 |
| app | `calibrate(y_meas, model_out, sensor, 'Bias', spline, 'Scale', "full"|"diag", 'Events', T) → Results (lsqvce)`; `forceBudget(Model, ep, r, v, q) → table`; `densityInversion(…)` | ○ | 19, 30, 32 |

Counting: 8 root + 10 io + 6 env + 7 geom + 8 kernel + 4 sensor + 6 qc + 3 app = **52 units** (`buildMex` added 2026-09-07), of which 25 ☑/interface, 27 ○ stubs (the 39th stub row, `Sensor`, became the interface; LP-2 turned seven io stubs into units on 2026-09-06). Every unit has one catalogue row, one file, one `Contents.m` line, and — when ☑ — seven test categories or a reasoned exemption.

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

### LP-4 — Geometry and attitude
WP-09 `boxWing`, cylinder/sphere kinds (the S-band boom) · WP-10 `selfShadow` raycast + LUT with a convergence sweep over `lutResolution` (T4) · WP-11 `AttitudeProvider` over `interpolateQuaternions`, `checkAttitudeSense` by physics (nadir, line of sight) · WP-28 `thermalState` steady-state and transient (oracle: [W] Chapter 4 figures re-digitised; energy balance metamorphic: total absorbed = total emitted in steady state).

### LP-5 — Kernel and quality control
WP-13 `gsiCoefficients` Sentman exact (oracle: Python kinetic-theory integration in `tools/pyxcheck/sentman_check.py`; the v1 ×8 defect as a regression test; flat-plate limits s→∞) · WP-26 CLL, Maxwell, Schamberg, `accommodation` SESAM (oracle Pilinski 2013 Fig. 3 values) · WP-12 `radiationPanelForce`, `thermalRecoil` (oracle: [V] monthly SRP/ERP means; [S] Table 3 RMS per axis within 20 %) · WP-24 `bodyTorques` (oracle: gravity-gradient closed form for a diagonal inertia; magnetic torque m×B) · WP-31 `qc`: `residuals`, `detectImpulses` (oracle: synthetic CRN-shaped pulses at known epochs, THR1B priors from the real file), `detectJumps` (oracle: planted steps, PELT vs brute-force segmentation on short series), `screen` (metamorphic: time shift, affine scaling, reversal). Accuracy: kernel partials vs complex step 1e-8; impulse epoch to ±1 sample; jump amplitude to 3σ.

### LP-5x — Cross-validation tooling (WP-33, WP-34; owner rulings C1(ii), C2(i), C3 of 2026-08-27)
**WP-33 `tools/pyxcheck/ngfkernel/`** — a numpy mirror of the FROZEN KERNEL ONLY: `gsiCoefficients`, `dragPanelForce`, `radiationPanelForce`, `thermalRecoil`, `assemble`. Written from the equations, not transcribed from the MATLAB, or it inherits the author's misunderstandings and stops being an oracle (T9). No mirror of io, env, geom, Model, sensor, qc. **Frozen when LP-5 is GREEN**: from then on a regression oracle, not a maintained second product. Rationale: LP-6b ports this same kernel to C, so both ports are then checked against one authority.

**WP-34 GROOPS bridge** — GROOPS is used as a BLACK BOX (GPLv3: no source is transcribed into this toolbox in any language). `SimulateAccelerometer` with `miscAccelerations:atmosphericDrag` (Sentman 1961 in the Moe & Moe 2005 analytical form, `considerTemperature`, `considerWind` HWM14) and `miscAccelerations:radiationPressure` (SRP + albedo/longwave grids + TRP) over the real 2020-01-01 GRACE-C day. Outputs are frozen as CSV under `tests/fixtures/groops/` and committed; neither the CI nor the owner's host needs a GROOPS installation. **Independence caveat:** GROOPS cites Wöske et al. (2019) for radiation and TRP — the same source as reference [W] — so for SRP/ERP/TRP this is an implementation comparison, not an independent oracle; the drag path (Moe & Moe analytical vs our exact Sentman) is a genuine cross-check.

Side effect: `gracefoReading.writeGroopsInstrument` / `readGroopsInstrument` document that their header conventions were never verified against a real GROOPS installation. WP-34 verifies them and reports back to that toolbox.

### LP-6 — Model and partials
WP-14 `Model.accel/torque` (integration requirement I2 from v0.1.1 kept: the kernel returns per-panel forces, never sums internally) · WP-15 `series` chunked/parallel with bit-identity test · WP-16 `partials` (r, v, Cd, Cr, α) · WP-17 `rhs()` and `budget` hook for gjToolbox · WP-30 `forceBudget`. Oracles: [S] and [V] magnitudes; ACT1B along-track drag at solar minimum ~5e-8 m/s².

### LP-6b — Single-epoch C/MEX
WP-17b transcribes the frozen kernel (GSI + radiation panel loop + assemble) to C, links `rsApiRotation`, parity suite vs MATLAB at 1e-14; only after LP-5/6 are GREEN.

### LP-7 — Sensors and applications
WP-29 `sensor.Electrostatic` (bias, 3×3 scale, CoM coupling from ω, ω̇ via `angularVelocityFromRotationRate`, CRN filter; oracle: [J] §3 and the GRACE-C CoM calibration of Hu 2022), `AtomInterferometer` (Δφ = k_eff a T², dead time, 2π ambiguity, rotation input; oracle Lévèque 2022 CARIOQA numbers), `DragFree` · WP-19 `calibrate` through lsqvce (oracle: recover the published GRACE-C along-track scale ≈1.05 and bias order from real ACT1B; recover planted parameters from synthetic data to their σ) · WP-32 `densityInversion`.

### LP-R — Release
WP-21 documentation build (`tools/buildDocs.m` regenerates `doc/html` from help blocks and the catalogue; `scanHelp.m` enforces the template; a static scan refuses Signal/Statistics toolbox calls) · WP-22 demos, benchmarks, `.mltbx` packaging, CI. No compatibility layer (F2).

---

## Part H — Obligations and open items

| # | Item | Owner | Re-entry |
|---|---|---|---|
| O-1 | `qc.hampel` loop is 382× slower than `movmedian` (0.64 s/day, acceptable); first candidate for the density/MEX session or a two-pass movmedian formulation with exact window-MAD | LP-5 | when a day-of-data screen exceeds 5 s end to end |
| O-2 | Handbook normals are printed to six digits (norm 1 − 5e-7); `MacroModel` tolerates 1e-5 and documents it; decide whether kernels renormalise (they currently do not) | LP-4 | before WP-13 |
| O-3 | CERES SYN1deg variable names (B7) | owner supplies a file | before WP-27 |
| O-4 | `gjToolbox` not on the path on the owner's host; `rhs()` contract to be confirmed against `gjToolbox.Dynamics` once L3 of that toolbox exists | owner | before WP-17 |
| O-5 | `referenceSystem` MEX reported "present but not loadable" by its own setup on 2026-08-27; **2026-09-06: setup reports `mex current` (gateway built 2026-08-05, SOFA 20231011, OpenMP 32 threads), so it loads on the host**; open only for the CI runner | owner | before LP-6b |
| O-6 | Exact TRP oracle: [S] Table 3 gives RMS per axis for GRACE-A 2008; a GRACE-FO 2020 oracle needs [V]'s 2025 update or ACT1B residual analysis | LP-5 | before WP-12 |
| O-7 | `Options.sun.source = "rs.de440"` requires furnished kernels; `radiationSources` must raise `depKernel` with the `downloadKernels` hint | LP-3 | WP-08 |
| O-8 | `readMAG1B`: no MAG1B file exists on the host (searched 2026-09-06), so there is no oracle and the unit stays a stub; the owner supplies a real MAG1B day (2020-01-01 GRACE-C) | owner | before WP-23/24 |
| O-9 | Dst and F30 are in none of the three space-weather sources in hand (CelesTrak SW-All; SET SOLFSMY, DTCFILE); DTM2020 (F30, Kp) and any Dst-driven model need another source | LP-3b | before WP-25 |
| O-10 | `referenceSystem.interpolateEop` / `rotationItrsToGcrs` refuse the shipped whole `eopc04_20.1962-now` (pre-1972 rows; identifier `leapSecondTable:beforeUtcEpoch`): report to the referenceSystem owner; NonGravForces documents the constraint (B.1, `readGNV1B` help) and lets the foreign error propagate | owner | before LP-R |
| O-13 | `srp.eclipse = "solaarsCF"` needs the SOLAARS-CF coefficient table (Robertson 2015 / Zhang et al. 2019), which is in no source on the host; `eclipseFactor` refuses it with `notImplemented` until the table is transcribed from the paper AND checked against a published curve | owner supplies the paper | before LP-6 (integrator C¹ requirement) |
| O-14 | MEX build on the CI runner: `buildMex` locates the MinGW-w64 support package itself when `MW_MINGW64_LOC` is unset (measured 2026-09-07: the host had the package at `%ProgramData%\MATLAB\SupportPackages\R2026a\3P.instrset\mingw_w64.instrset` and MATLAB did not know); `ciGate` builds before the suite; a runner without the package turns the gate red by design | owner | first CI run of LP-3 |
| O-15 | `radiationSources` with `sun.source = "rs.de440"` is delivered as code (mice branch through `referenceSystem.geocentricBodyState`) and verified only for its `depKernel` refusal; the furnished-kernel path is ◐ until a kernel is available to the suite on every tier | LP-3b | with O-7 |
| O-12 | `gracefoReading.readL1bMulti` 1.0.0 raises `MATLAB:badsubscript` (line 156, `time(1,1)` on an empty concatenation) when TimeRange excludes every record, while `readL1b` returns clean empties — report to the gracefoReading owner; until fixed, `io/private/l1bLoad` reads one file through `readL1b` and only several through `readL1bMulti`, so `dataEmpty` is guaranteed for one file and documented as unreachable for several | owner | before LP-R |
| O-11 | `doc/html/FunctionList.html` is regenerated by `tools/buildFunctionList.py` and checked for staleness by `auditPackage.py`; the five other pages carry hand-edited version footers (checked for the version string only); the full documentation build is WP-21 (LP-R) | LP-R | WP-21 |

---

## Part I — Decision log

| Date | Decision | Reason |
|---|---|---|
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
