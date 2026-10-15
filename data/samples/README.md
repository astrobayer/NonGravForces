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
