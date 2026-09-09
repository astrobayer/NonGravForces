# CLAUDE.md — working notes for NonGravForces

Orientation file for an AI assistant working in this repository.
`HANDOVER.md` is the single source of truth; this file is the short form.
Read `HANDOVER.md` Part A, then the level prompt in Part G for the level being worked.

Created 2026-09-06 by Claude Opus 5 (Anthropic); updated 2026-09-08 (LP-7 opened, over the same day's LP-6b) by Claude Opus 5. Toolbox version 0.8.0.

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

**The rule that governs every other rule:**
> NonGravForces computes forces and torques and nothing else.

No gravity field, no equations of motion, no integrator, no least-squares code.

---

## 2. Architecture

Nine layers; a layer calls only lower layers plus the four dependencies.

```
L8 app      calibrate (lsqvce bridge), forceBudget, densityInversion
L7 qc       residuals, hampel, detectImpulses, detectJumps, eventTable, screen
L6 sensor   Sensor (abstract), Electrostatic, AtomInterferometer, DragFree
L5 Model    accel | torque | series | partials | rhs | budget   <- single point of access
L4 kernel   gsiCoefficients, accommodation, dragPanelForce, radiationPanelForce,
            thermalRecoil, assemble, bodyTorques, partials, singleEpoch (C/MEX)
L3 geom     MacroModel, graceFO, boxWing, selfShadow, thermalState,
            AttitudeProvider, checkAttitudeSense
L2 env      density, wind, eclipseFactor, earthRadiance, magneticField, radiationSources
L1 io       readACC1B/SCA1B/GNV1B/THR1B/MAS1B/MAG1B, massModel, spaceWeather,
            importCERES, igrfCoefficients
L0 core     const, toolboxVersion, dependencyReport, fail, provenance, Options, buildMex
```

53 public units, catalogued in `doc/functionlist.csv` (the authority; CI checks
both directions against the files on disk, and this sentence's number against it).

### Dependencies — called by qualified name, never copied, wrapped or aliased

| Toolbox | Owns |
|---|---|
| `referenceSystem` ≥ 2.1 (2.1.216) | time scales, EOP, ITRS↔GCRS with derivatives, geodesy, quaternions, attitude interpolation, Sun/Moon states |
| `gracefoReading` ≥ 1.0.1 (1.0.1) | L1B ASCII reading, quality flags, GROOPS instrument export/import. The PATCH minimum is deliberate and is the one exception to VERSIONING.md's MAJOR.MINOR rule (HANDOVER O-12, O-25) |
| `gjToolbox` (odeSuite) | Gauss-Jackson integration. **At level L0: no integrator and no dynamics interface yet, and `Model.rhs()`'s 4-argument handle does not match the 3-argument one its validator documents — HANDOVER B.3, O-4** |
| `lsqvce` ≥ 5.3 (5.3) | batch least squares, VCE, ARMA whitening; `app.calibrate` delegates the adjustment through `fit`/`layout`/`Options` (O-24) |

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
  `tNgfMex/singleEpochEndToEndPrecision` is where that is checked.
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

## 5. Current status (2026-09-07, after LP-5x and the debt sweep)

- **276/276 GREEN on the CI runner** (2026-09-08, LP-7 WP-30, R2026a Update 5, 219.6 s; predicted
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
- **4 documented stubs** raising `NonGravForces:<unit>:notImplemented`.
- LP-6 added the six `Model` evaluation methods (`accel`, `torque`, `series`, `partials`,
  `rhs`, `budget`), `env.wind` (co-rotation) and `kernel.partials` (analytic derivatives).
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

**Done 2026-09-07: LP-5 (kernel + qc)** except `partials` (LP-6), SESAM (O-16),
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
The Coriolis phase at orbital rates is **6.7e7 rad**, an apparent 16.7 m/s^2: not a correction to
the measurement but larger than everything the instrument exists to see.

`sensor.DragFree` (WP-29) closes the level: the residual a loop FAILS to remove, as a cascade of
bilinear high-pass sections with actuator saturation and a masked start-up transient.

**LP-7 IS COMPLETE.** The four remaining stubs are blocked on data the project does not hold, so
they are not level debt. **Next is LP-R (release).** Two owner rulings of 2026-09-09 shape it: the
toolbox is **MIT** (`LICENSE`), and **O-29 is accepted** — 88 ms per epoch in the standard preset
ships as a stated cost, not a blocker. **O-34 is the one thing that must be settled before the
repository goes public**: the vendored NRLMSISE-00 source carries no redistribution grant, and
`THIRD-PARTY-NOTICES.md` recommends removing it in favour of a user-supplied copy.

**LP-5x cross-validation** (owner rulings of 2026-08-27):
- WP-33: standard-library mirror of the **kernel layer only** under
  `tools/pyxcheck/ngfkernel/`, written from the equations (not transcribed from
  the MATLAB), **frozen once LP-5 is GREEN** as a regression oracle.
- WP-34: GROOPS as a **black box** (GPLv3 — no source transcribed in any
  language). `SimulateAccelerometer` output for the real day is frozen as CSV
  under `tests/fixtures/groops/` and committed; neither CI nor the host needs a
  GROOPS installation.

**Obligations**: `HANDOVER.md` Part H holds them, currently O-1 to O-31 with an owner
and a re-entry condition each. This file deliberately does NOT list them: the copy that
used to live here had drifted to O-15 of 21 and showed two closed items as open, and no
check compared the two (R1). Read Part H.

Housekeeping: the repository has no LICENSE file, and the vendored NRLMSISE-00 sources
carry no redistribution grant (O-22, the owner's to settle before the repository is public).

**Version rule (decided 2026-09-06):** the numeric fields of `toolboxVersion`
are the authority and `semantic` is composed from them; minor = levels
delivered, patch = HANDOVER amendment. `auditPackage.py` compares every copy
(Contents.m, README, doc/html footers) and checks `FunctionList.html` against
`tools/buildFunctionList.py`; run both before every push.

---

## 7. Known pitfalls

Every one of these was an actual defect found in a predecessor or a dependency.

- **MAS1B is cold-gas tank mass** (~15 kg per tank, ~30.5 kg total), **not
  spacecraft mass**. Use the handbook constant: 601.214 kg (GRACE-C),
  601.209 kg (GRACE-D). Misuse gives 20–40× accelerations. `io.massModel`
  defaults to the constant and asserts > 500 kg.
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
  projections. The theory document §2.3 mislabels them (to be corrected at LP-R). Three
  independent derivations reached the same conclusion in LP-5x. The limit that settles it:
  edge-on, Cp = 1/(2s^2) and Ctau = l/(s*sqrt(pi)).
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
- **Never claim bitwise across a change in BATCH SIZE.** Vectorised products and norms associate
  by batch size, so K = 1 and K = 10 are not the same arithmetic even when they are the same
  formula. Met three times now: `pagemtimes` in `Model.evaluate`, the chunked RMS in
  `forceBudget`, and one epoch vs an arc in `densityInversion`. Bitwise IS claimable across a
  change of grouping that leaves the arithmetic alone — two models stacked, a maximum of maxima.
- **A numerical bound written from an estimate is a red gate waiting to happen.** Three sessions
  in a row have lost a run to one. Where cancellation governs, COMPUTE the bound in the test from
  `eps` of the terms involved instead of writing a literal. Where the model is exactly linear, a
  LARGER finite-difference step is strictly better, the opposite of the usual rule.
- **`skew(omega) = R' Rdot` is the sensor's convention**, and
  `referenceSystem.angularVelocityFromRotationRate(R, Rdot, 'Resolution', "reference")` returns its
  NEGATIVE (ratio -1.0, measured). Reading the skew entries off `R' Rdot` without projecting the
  symmetric part out leaves the differencing error in (6e-7 relative).
- **Never assert a host-dependent absolute** (rule T5). A 500 us per-epoch bound written from a
  guess reddened the gate on a runner that measured 901.8 us, in the same session that quoted the
  rule. Assert the ratio; record the absolute with a ceiling an order of magnitude above it.
- **When the MCP bridge is down, CI is still a MATLAB tier**: push the branch, let the `matlab`
  job run `ciGate`, read `tests/last_run.txt` from the run artefact. That is how LP-6b was verified.
