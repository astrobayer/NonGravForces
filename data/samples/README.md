# data/samples — real GRACE-FO Level-1B products (test oracles)

Six RL04 ASCII products of **2020-01-01, GRACE-C** from PO.DAAC
(https://podaac.jpl.nasa.gov/GRACE-FO), taken over unchanged from the owner's
predecessor toolbox (`NonGravToolbox/data/samples`, 2026-08-01):

| File | Records | Content |
|---|---|---|
| `ACT1B_2020-01-01_C_04.txt` | 600 (first 10 min) | calibrated linear + angular accelerations, SRF |
| `GNV1B_2020-01-01_C_04.txt` | 600 (first 10 min) | position and velocity, Earth-fixed (`coord_ref = E`) |
| `SCA1B_2020-01-01_C_04.txt` | 600 (first 10 min) | star-camera quaternions ICRF→SRF, scalar first |
| `THR1B_2020-01-01_C_04.txt` | 471 (complete day) | thruster activations (on-times per thruster, ms) |
| `MAS1B_2020-01-01_C_04.txt` | 24 (complete day) | cold-gas TANK masses (~15 kg per tank) — NOT spacecraft mass |
| `MAG1B_2020-01-01_C_04.txt` | 600 (first 5 min, 0.5 s) | measured field in SRF (µT) and six magnetorquer branch currents (mA) |

The three 1 Hz products and MAG1B (0.5 s) are the first 600 records of the day
with `num_records` rewritten to 600; every data line is as published. MAG1B was
cut from the GFZ ISDC archive on 2026-09-09 by the same recipe. MAS1B and THR1B are
byte-identical to `gracefoReading/tests/realdata`.

These files are the oracle of `tests/tNgfIo.m` (HANDOVER Part I, 2026-09-06).
Replacing an excerpt with the full-day file needs no code change: readers and
tests read record counts from the file. Tests never skip when a file is
missing; they fail, because a skipped test is not evidence.

## multiday/ — three consecutive days for the multi-arc calibration (N-09)

`ACT1B`, `GNV1B` and `SCA1B` of **2020-01-01, -02 and -03, GRACE-C**, RL04 ASCII, each the
**first 1800 records (30 min)** of the day with `num_records` rewritten to 1800; every data line
is as published. Cut on 2026-10-15 from the daily archives of GFZ ISDC
(`https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04/2020/`), anonymously available:

| Archive | Bytes | SHA-256 |
|---|---|---|
| `gracefo_1B_2020-01-01_RL04.ascii.noLRI.tgz` | 146773459 | `f14fec9ab9e30ce3cba95b14c04a5e03fc5ab66815076a17b5e72eba0ef6557e` |
| `gracefo_1B_2020-01-02_RL04.ascii.noLRI.tgz` | 147026710 | `b13e61ea46a2baab3cf97023ad4f8a86f5e7677faec21b2e68f0631206f2ece7` |
| `gracefo_1B_2020-01-03_RL04.ascii.noLRI.tgz` | 147500819 | `116453c6ef4ff09ec3f98df5276c7fd298438a88365ad2ff84dd476013203a87` |

The 2020-01-01 excerpts here are longer than the 600-record ones above, which stay as they are:
`tests/tNgfIo.m` reads those. These are the oracle of `tNgfCalibrate/calibrateCampaignRealDays`.

Added for N-11a (2026-10-15), from the same 2020-01-01 archive: `AHK1B_2020-01-01_C_04.txt`, the first
**30 min (18002 records)** of accelerometer housekeeping, variable columns per record as `prod_flag` sets
them; and `MAG1B_2020-01-01_C_04.txt`, the first **3600 records (30 min at 0.5 s)**; both with `num_records`
rewritten and every data line as published. They are the oracle of `tNgfCalibrate/calibrateRegressorsRealDay`.
On that day the sensor-unit temperatures TSU span only 0.022 C and TICUR reads -273.2 C (a dead channel).

## itsg/ — TU Graz operational orbit products (N-12 attitude; an external force series for O-19)

Excerpts of the Institute of Geodesy's operational products, `https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational/<satellite>/<product>/<year>/`,
fetched 2026-10-26 by `.github/workflows/fetch-data.yml` (run 40361879447; `tools/fetchData.py`). Free of charge for
any use; publications acknowledge TU Graz (their `README.txt`). GRACE-1/-2 are GRACE-A/-B, GRACEFO-1/-2 are GRACE-C/-D.

Each file is the **first 600 epochs (10 s sampling, 00:00:00 to 01:39:50 GPS)** of the day, GROOPS instrument format
(`groops instrument version=20200123`), with the epoch-count line after the `# ===` rule rewritten from 8641 to 600;
every data line is as published. Times are GPS, as MJD. Frames, from TU Graz's `README.txt`:

| Product | Content | Frame |
|---|---|---|
| `attitude` | quaternion, scalar first | rotation satellite reference frame -> celestial reference frame |
| `nonConservativeForces` | acceleration due to non-conservative forces, m/s^2 | satellite reference frame |
| `reducedDynamicOrbit` | position, velocity (and three zero acceleration columns), m, m/s | celestial reference frame |

TU Graz's README does not say whether `nonConservativeForces` is the calibrated accelerometer or a model; nothing in
this repository may assume either until that is settled with them. The GRACE-FO orbits are in their own CoM-realised
quasi-inertial frame (`README_referenceFrames.txt`), not referenceSystem's GCRS: compare after a frame check.

| Downloaded file (gz) | Bytes | SHA-256 |
|---|---|---|
| `GRACE-1_attitude_2008-01-01.txt.gz` | 438416 | `b6689a4f2c0fc1f28a5db9eba82dd28a2f4a8f6bf1a77d53b2b38f264446b386` |
| `GRACE-2_attitude_2008-01-01.txt.gz` | 437714 | `11f8319ce529933c3945fe75ad2027f24a5e77e6b4677c27b65904af01dcb1fe` |
| `GRACEFO-1_attitude_2020-01-01.txt.gz` | 439031 | `6c845b358f91175aad28f87f70ef4605e0b0abe253b015fc88fc24c5dcd5ba72` |
| `GRACEFO-1_attitude_2020-01-02.txt.gz` | 438820 | `c7614304a1b30ef15c7ce2660f9c3c6a7690c0bbf3ac91cfc914b06b803d2f89` |
| `GRACEFO-1_attitude_2020-01-03.txt.gz` | 438903 | `8776843d45c936a32cd8fc461a63e7e8371cc8e6c785ef400add4cc7a7582d4e` |
| `GRACEFO-2_attitude_2020-01-01.txt.gz` | 437144 | `523f83d00bff30dc8193ed06e876cb24032e14da4b8d32c00913a3480ce2e1dd` |
| `GRACE-1_nonConservativeForces_2008-01-01.txt.gz` | 341164 | `4e720c8b730435c94f17a334e8ec291e414673f557f2fd6b2acdcb89e295b985` |
| `GRACE-2_nonConservativeForces_2008-01-01.txt.gz` | 340654 | `8326d5563d5031f725f5bb858084df9b8b4266dc6f7372a5c9db79994c83849b` |
| `GRACEFO-1_nonConservativeForces_2020-01-01.txt.gz` | 340331 | `57023e3f6383e843a325b597b81c5e24624397fc4e9b135b481f90721338034b` |
| `GRACEFO-1_nonConservativeForces_2020-01-02.txt.gz` | 340232 | `5528e19a783b924add176a849d2904217437d60c7e7f517efff4ddda0b2b4914` |
| `GRACEFO-1_nonConservativeForces_2020-01-03.txt.gz` | 340753 | `3a8316726b8fc71de60492affe4af63e4789a15e4f451e279e59b4595a4430c9` |
| `GRACEFO-2_nonConservativeForces_2020-01-01.txt.gz` | 339753 | `3d3490e1ed629867a6481cd8c80ae574c3a0dd6f36f5625ecd4a76dba9b95723` |
| `GRACE-1_reducedDynamicOrbit_2008-01-01.txt.gz` | 661333 | `40813ac1aeb1ad2177e1b19e48ced675d272b240ba88247a34289d8e4ca4801b` |
| `GRACE-2_reducedDynamicOrbit_2008-01-01.txt.gz` | 661454 | `e3d9d7463f5e7ecfa34acae3bc09f64619151d00971242883d2b64d448a51aa6` |
| `GRACEFO-1_reducedDynamicOrbit_2020-01-01.txt.gz` | 659073 | `32608c82721ac7ca8e97291687475a25f880738e7009c3319fe916e425868f7b` |
| `GRACEFO-1_reducedDynamicOrbit_2020-01-02.txt.gz` | 658882 | `24d9065d9895388a41a4dc68978504d6b2e64e7bb8a826578822181594a1daf1` |
| `GRACEFO-1_reducedDynamicOrbit_2020-01-03.txt.gz` | 658883 | `dc46ecf50d9841629a12a9ea5c47e00dbf6a68e86bb2448e5ec193260072524c` |
| `GRACEFO-2_reducedDynamicOrbit_2020-01-01.txt.gz` | 659035 | `e90356d8869d8662dbaf8832b3ef57bd7139a4a69aeaf04e7c20999f449426d4` |
