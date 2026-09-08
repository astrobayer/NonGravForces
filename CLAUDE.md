# CLAUDE.md — working notes for NonGravForces

Orientation file for an AI assistant working in this repository.
`HANDOVER.md` is the single source of truth; this file is the short form.
Read `HANDOVER.md` Part A, then the level prompt in Part G for the level being worked.

Created 2026-09-06 by Claude Opus 5 (Anthropic); updated 2026-09-08 (O-12 discharge) by Claude Opus 5. Toolbox version 0.5.3.

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
            thermalRecoil, assemble, bodyTorques, partials
L3 geom     MacroModel, graceFO, boxWing, selfShadow, thermalState,
            AttitudeProvider, checkAttitudeSense
L2 env      density, wind, eclipseFactor, earthRadiance, magneticField, radiationSources
L1 io       readACC1B/SCA1B/GNV1B/THR1B/MAS1B/MAG1B, massModel, spaceWeather,
            importCERES, igrfCoefficients
L0 core     const, toolboxVersion, dependencyReport, fail, provenance, Options, buildMex
```

52 public units, catalogued in `doc/functionlist.csv` (the authority; CI checks
both directions against the files on disk, and this sentence's number against it).

### Dependencies — called by qualified name, never copied, wrapped or aliased

| Toolbox | Owns |
|---|---|
| `referenceSystem` ≥ 2.1 (2.1.216) | time scales, EOP, ITRS↔GCRS with derivatives, geodesy, quaternions, attitude interpolation, Sun/Moon states |
| `gracefoReading` ≥ 1.0.1 (1.0.1) | L1B ASCII reading, quality flags, GROOPS instrument export/import. The PATCH minimum is deliberate and is the one exception to VERSIONING.md's MAJOR.MINOR rule (HANDOVER O-12, O-25) |
| `gjToolbox` (odeSuite) | Gauss-Jackson integration; `Model.rhs()` hands it an acceleration handle |
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
- **MATLAB is normative.** A C/MEX port (LP-6b, hand-written, no MATLAB Coder)
  and the standard-library kernel mirror (WP-33) are verified *against* MATLAB, never
  the reverse.

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

- **244/244 tests GREEN** (2026-09-07, R2026a Update 4, Tier A host run; predicted 243,
  244 ran — reconciled in HANDOVER Part I), Code Analyzer 0 findings, gate TRUE;
  log in `tests/last_run.txt`.
- **245/245 tests GREEN** (2026-09-08, self-hosted Windows runner, `main` after PR #16;
  **predicted 245, 245 ran** — exact, 244 plus the one robustness method
  `tNgfIo/readACC1BMultiEmptyRobustness`), 0 audit shortfalls, gate TRUE. This is the
  0.5.3 (O-12 discharge) run report, so 0.5.3 is ☑, not ◐: it was written on Tier B and
  verified on the runner afterwards.
  LP-5x was 238, LP-5 231, LP-4 169, LP-3 142, LP-2 116 (2026-09-06), LP-0 74 (2026-08-27).
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
- **12 documented stubs** raising `NonGravForces:<unit>:notImplemented`.
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

**Next: LP-6 (Model and partials)** — WP-14 `Model.accel/torque`, WP-15 `series`
(chunked, parfor-optional, bit-identical), WP-16 `partials` and `kernel.partials`
(analytic, verified by complex step), WP-17 `rhs()` for gjToolbox; oracles in
HANDOVER Part G. Then LP-4 (geom), LP-5 (kernel + qc), LP-5x, LP-6 (Model/partials),
LP-6b (C/MEX), LP-R (release).

**LP-5x cross-validation** (owner rulings of 2026-08-27):
- WP-33: standard-library mirror of the **kernel layer only** under
  `tools/pyxcheck/ngfkernel/`, written from the equations (not transcribed from
  the MATLAB), **frozen once LP-5 is GREEN** as a regression oracle.
- WP-34: GROOPS as a **black box** (GPLv3 — no source transcribed in any
  language). `SimulateAccelerometer` output for the real day is frozen as CSV
  under `tests/fixtures/groops/` and committed; neither CI nor the host needs a
  GROOPS installation.

**Obligations**: `HANDOVER.md` Part H holds them, currently O-1 to O-24 with an owner
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
