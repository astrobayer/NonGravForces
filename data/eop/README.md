# data/eop — Earth orientation excerpt (test oracle)

`finals2000A_2019-09_2020-03.all`: the IERS Bulletin A `finals2000A.all` rows
for 2019-09-01 .. 2020-03-31 (213 days, all flagged I = final), cut by date
from the owner's copy (referenceSystem data folder, fetched 2026-08). Read with
`referenceSystem.readEop`; it rotates the 2020-01-01 GNV1B sample bitwise
identically to the full series (measured 2026-09-06).

Why an excerpt is committed: referenceSystem refuses an EOP series with rows
before 1972-01-01 (HANDOVER O-10), so the natural choice, the whole C04
1962-now file, cannot be handed to `readGNV1B`; and a test must never depend on
the host's data folder. referenceSystem ships the full `finals2000A.all` in its
own `data/` for production use.
