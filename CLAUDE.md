# CLAUDE.md — working notes for NonGravForces

Orientation file for an AI assistant working in this repository.
`HANDOVER.md` is the single source of truth; this file is the short form.
Read `HANDOVER.md` Part A, then the level prompt in Part G for the level being worked.

Created 2026-09-06 by Claude Opus 5 (Anthropic). Toolbox version 0.1.0.

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
L0 core     const, toolboxVersion, dependencyReport, fail, provenance, Options
```

51 public units, catalogued in `doc/functionlist.csv` (the authority; CI checks
both directions against the files on disk).

### Dependencies — called by qualified name, never copied, wrapped or aliased

| Toolbox | Owns |
|---|---|
| `referenceSystem` ≥ 2.1 (2.1.214) | time scales, EOP, ITRS↔GCRS with derivatives, geodesy, quaternions, attitude interpolation, Sun/Moon states |
| `gracefoReading` ≥ 2.1 (2.3.0) | L1B ASCII reading, quality flags, GROOPS instrument export/import |
| `gjToolbox` (odeSuite) | Gauss-Jackson integration; `Model.rhs()` hands it an acceleration handle |
| `lsqvce` 1.0 | batch least squares, VCE, ARMA whitening; `app.calibrate` delegates the adjustment |

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
  and the numpy kernel mirror (WP-33) are verified *against* MATLAB, never the
  reverse.

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
`main` is protected: PR required, both checks required, strict, no force push.
Work goes on a branch → PR → green → merge → `git pull` on the host.

---

## 5. Current status (2026-09-06)

- **74/74 tests GREEN**, Code Analyzer 0 findings, gate TRUE
  (last verified run 2026-08-27, R2026a Update 4; log in `tests/last_run.txt`).
- Implemented (☑): L0 core, `geom.MacroModel`, `geom.graceFO`, `qc.hampel`,
  `qc.eventTable`, `sensor.Sensor` interface, test infrastructure, CI.
- **39 documented stubs** raising `NonGravForces:<unit>:notImplemented`.
- Repository `github.com/astrobayer/NonGravForces` (the account was renamed
  from `matweigelt` on 2026-09-06; `referenceSystem` and `gracefoReading` moved
  with it), local clone at
  `C:\Users\matth\Documents\MATLAB\NonGravForces` with the dependencies as
  siblings.

---

## 6. Open tasks

**Next: LP-2 (Data I/O)** — WP-18 readers over `gracefoReading`, WP-04
`spaceWeather`, `massModel`. Blocked only on real L1B files (ACT1B, SCA1B,
GNV1B, MAG1B for 2020-01-01 GRACE-C; MAS1B and THR1B already exist in
`gracefoReading/tests/realdata`). Oracles: GNV1B round trip < 1 mm, exact
time-tag round trip, Greenwich-angle test of the rotation direction.

Then LP-3 (env), LP-4 (geom), LP-5 (kernel + qc), LP-5x, LP-6 (Model/partials),
LP-6b (C/MEX), LP-R (release).

**LP-5x cross-validation** (owner rulings of 2026-08-27):
- WP-33: numpy mirror of the **kernel layer only** under
  `tools/pyxcheck/ngfkernel/`, written from the equations (not transcribed from
  the MATLAB), **frozen once LP-5 is GREEN** as a regression oracle.
- WP-34: GROOPS as a **black box** (GPLv3 — no source transcribed in any
  language). `SimulateAccelerometer` output for the real day is frozen as CSV
  under `tests/fixtures/groops/` and committed; neither CI nor the host needs a
  GROOPS installation.

**Obligations** (`HANDOVER.md` Part H): O-1 `hampel` speed · O-2 handbook normal
renormalisation · O-3 CERES SYN1deg variable names · O-4 `gjToolbox` not on the
host path · O-5 `referenceSystem` MEX "present but not loadable" (needed for
LP-6b) · O-6 exact TRP oracle for GRACE-FO 2020 · O-7 DE440 kernels for
`Options.sun.source = "rs.de440"`.

Housekeeping: the repository has no LICENSE file yet.

**Open defect to confirm:** `toolboxVersion` now returns `semantic = "0.1.0"`
while the numeric fields still read `major = 2, minor = 0, patch = 0`. One fact,
two authorities (violates R1) — decide which is normative and make the other
follow.

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
- **Sentman**: the shear term divides by s², not s (v1 was ×8 wrong at LEO
  s≈8); the ambient gas is atomic O, not air — `R = 287 J/kg/K` is wrong by
  √(29/16) = 1.35 at 500 km.
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
