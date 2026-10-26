# data/spaceweather — solar and geomagnetic index excerpts (test oracles)

Excerpts for **2019-09-01 .. 2020-03-31** (213 days) cut by date from the
owner's copies (`NonGravToolbox/data/spaceweather`, downloaded 2026-08); header
lines kept, data lines unchanged:

| File | Source | Serves |
|---|---|---|
| `SW-All_2019-09_2020-03.csv` | CelesTrak, https://celestrak.org/SpaceData/SW-All.csv | family `celestrak`: F10.7 (observed, adjusted, 81-day means), 3-hourly Kp and ap, daily Ap |
| `SOLFSMY_2019-09_2020-03.TXT` | Space Environment Technologies, release 8_1_0 (2026-07-31), https://sol.spacenvironment.net/JB2008/indices/SOLFSMY.TXT | family `jb2008`: F10, S10, M10, Y10 and their 81-day centred means |
| `DTCFILE_2019-09_2020-03.TXT` | Space Environment Technologies, https://sol.spacenvironment.net/JB2008/indices/DTCFILE.TXT | family `jb2008`: hourly exospheric temperature change DTC |

`NonGravForces.io.spaceWeather` identifies the family from the file CONTENT
(never from the name). The 81-day means are the published columns, so an
excerpt carries the same values as the full file. Dst and F30 are in none of
these sources (HANDOVER O-9); they came separately, below.

## Dst and F30 (O-9), fetched 2026-10-26

Fetched by `.github/workflows/fetch-data.yml` (`tools/fetchData.py`, run 40361879447) on a GitHub-hosted runner,
because the development container cannot reach these hosts. **No reader exists for either yet**:
`io.spaceWeather` identifies the three families above only, and these files are the oracle a reader must be
written against.

| File | Source | Content |
|---|---|---|
| `dst0712.for.request`, `dst0801.for.request`, `dst1912.for.request`, `dst2001.for.request` | WDC for Geomagnetism, Kyoto, **final** Dst, `https://wdc.kugi.kyoto-u.ac.jp/dst_final/<yyyymm>/dst<yymm>.for.request` | hourly Dst (nT), WDC exchange format, one record per day (24 hourly values and the daily mean); complete months, unchanged |
| `radio_flux_adjusted_observation_2007-09_2008-03_2019-09_2020-03.txt` | CLS, `ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation/radio_flux_adjusted_observation.txt` | daily f30, f15, f10.7, f8, f3.2 (sfu) adjusted to 1 AU, each with flare-corrected value, precision and processing flag (1 provisional, 2 interpolated, 4 flare-corrected) |
| `radio_flux_absolute_observation_2007-09_2008-03_2019-09_2020-03.txt` | CLS, same folder, `radio_flux_absolute_observation.txt` | the same, not adjusted to 1 AU |

The CLS files keep their header unchanged and the data lines of 2007-09-01..2008-03-31 and 2019-09-01..2020-03-31
(213 + 213 daily lines, every line as published). CLS fills gaps by interpolation and flags it: on
the GRACE-FO fixture days 2020-01-01..03 f30 is flagged 2 (interpolated), because Nobeyama's own daily file
(`TYKW-NoRP_dailyflux.txt`) has no 1 GHz value on those days. DTM2020 is driven by the CLS series.

| Downloaded file | Bytes | SHA-256 |
|---|---|---|
| `dst0712.for.request` | 3752 | `7b3ff06cf0cc4560497e8ca9030d9b67f6f118e396ff40657daef86b0b7fc432` |
| `dst0801.for.request` | 3752 | `a09c56d390e89bad6d7dd89a3db7858503e52649f6087fd57aee594b8707b23b` |
| `dst1912.for.request` | 3752 | `8813b818ffaa2a04521be4ed0f225a8164e490c520b4d4c408dc1c0211f8e99a` |
| `dst2001.for.request` | 3752 | `f7545be24ea25aecb985bb5b73f760cd620d2897bbe041bcd6145761c45da483` |
| `radio_flux_adjusted_observation.txt` | 5260147 | `6423391a3bd0b9fb073170c0a896ff9b6c22f1db0e5ed6c554e7a3a358f5e422` |
| `radio_flux_absolute_observation.txt` | 5260151 | `50c8df9ab77ebe0be6387a0d596e21956e71f5016458d4252790d8d6514cf46f` |
