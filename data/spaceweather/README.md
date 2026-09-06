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
these sources (HANDOVER O-9).
