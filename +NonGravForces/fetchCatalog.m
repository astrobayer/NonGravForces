function [T, plan] = fetchCatalog(options)
%FETCHCATALOG Everything setupNonGravForces can fetch: spans, sizes, licence classes, and what a selection implies.
%
%   T = NonGravForces.fetchCatalog
%   [T, plan] = NonGravForces.fetchCatalog('PoolDays', "2020-01-01..2020-01-31", 'Items', ["itsg" "dst"])
%   [T, plan] = NonGravForces.fetchCatalog('PoolMonths', "all")
%
%   One row per selectable SOURCE (an id), not per file: the sources a caller can choose with
%   setupNonGravForces('Items', ...). The table says where each source goes, which setup option
%   governs it, the span the provider publishes, what one file costs in download and on disk,
%   and the licence class this repository has RECORDED for it. It reaches no network and reads
%   no file outside the toolbox; it is the catalogue an installer shows before anything is
%   fetched.
%
%   THE ONE AUTHORITY FOR THREE FACTS (R1). setupNonGravForces owns the URLs, pins and the fetch
%   itself; this unit owns (1) the ids, (2) the spans of the dated products, and (3) how a
%   selection is expanded. setup calls this unit with its own options and fetches exactly the
%   expanded plan, so "all" and "a..b" mean the same thing here and there.
%   tNgfFetchCatalog/fetchCatalogMetamorphic compares the counts of this table with the records
%   setup produces for the same selection, per id, so the two cannot drift apart unnoticed.
%
%   SELECTION SYNTAX (PoolDays, PoolMonths, SwarmMonths, DensityMonths). Each element is one of
%     "yyyy-mm-dd" or "yyyy-mm"                    a single day or month
%     "yyyy-mm-dd..yyyy-mm-dd", "yyyy-mm..yyyy-mm" an inclusive range
%     "all"                                        every day or month the provider publishes
%   Elements are united, sorted and de-duplicated. "all" for PoolDays is the union of TU Graz's
%   operational spans (GRACE, then GRACE-FO), about 8800 days and 25 GB.
%
%   THE TOLEOS GAP (DensityMonths). No TOLEOS mission flew 2017-11..2018-04 (GRACE-A's last day
%   2017-10-31, GRACE-FO 1's first 2018-05-29). "all" is the union of the three missions' own
%   spans: 2000-07..2017-10 and 2018-05..2026-04, 304 months, the gap months absent. An explicit
%   range that crosses the gap keeps its flown months and SKIPS the gap's, which hold no file; the
%   plan lists exactly the months kept, so nothing is dropped unseen ("2017-10..2018-05" is two
%   months). An element that selects no flown month is refused: a single gap month ("2017-11")
%   or a range inside the gap. A range reaching past the last TOLEOS month is refused too: that
%   is outside the published span, not a gap. Until 2026-11-09 "all" was expanded as the one
%   envelope 2000-07..2026-04 and refused itself at 2017-11.
%
%   LICENCE CLASSES. The class is what THIRD-PARTY-NOTICES.md records, and nothing more:
%     "open"        reuse allowed, commercial included, with attribution (CC BY 4.0, ESA's open data)
%     "academic"    research, academic and non-profit use only; commercial use prohibited
%     "unstated"    the package carries no licence and no grant; nothing may be inferred from that
%     "unreviewed"  provider terms apply and this repository has not reviewed them
%   An installer should ask for an acknowledgement before fetching "academic" or "unstated".
%   Column notice names the heading of THIRD-PARTY-NOTICES.md that records the class ("" for
%   "unreviewed"); tNgfFetchCatalog/fetchCatalogContract checks that every named heading exists.
%
%   WHERE EACH SOURCE GOES (column destination; owner, 2026-11-08: "All data go always to the
%   data pool. Only model code or libraries go to the toolbox."). Every DATA source is "pool": it
%   is written below the root NonGravForces.dataPool resolves ('DataPool', else DATAPOOL, else
%   READING_DATA_POOL), in the product folder every astrobayer toolbox shares (SPICE, EOP,
%   SpaceWeather, IGRF, CERES/EBAF, EarthRadiationPressure, GRACEFO/L1B/RL04, ITSG, Swarm, TOLEOS).
%   Only the four MODEL SOURCES (nrlmsise00, nrlmsis21, hwm14, dtm2020, with the coefficient files
%   they ship with) are "toolbox": they are compiled into its MEX. Until 2026-11-08 eop,
%   spaceweather, igrf, l1b-reference, ceres and de440 went into the toolbox's data/ folder.
%
%   SPANS AND SIZES ARE MEASURED, AND SAY SO. Column sizeBasis says where each figure comes from.
%   Measured on 2026-11-06 (the owner's host and pool, MATLAB webread): TU Graz's operational
%   products from 2002-04-04 to 2017-06-29 (GRACE) and from 2018-06-01 (GRACE-FO, last day
%   2026-09-30 on that date); one ITSG file 0.34-0.66 MB; a TOLEOS GRACE-FO month 39.9 MB; the
%   five GRACE-FO Level-1B members of one day for C and D 299 MB on disk. TOLEOS spans are the HAPI
%   headers of the committed excerpts (data/samples/toleos), Swarm's first day the header of
%   data/samples/swarm (2014-01-01). An "ongoing" product ends, for "all", on the last day of the
%   month before 'Today': later files may not be published yet, and setup reports each missing one.
%
% INPUT ARGUMENTS
%   Name   Size  Type  Units  Description
%   ----   ----  ----  -----  -----------
%   (none)
%
% NAME-VALUE ARGUMENTS
%   Name           Size  Type      Units  Description
%   ----           ----  ----      -----  -----------
%   PoolDays       1xN   string    -      Reference days for FetchPool (Dst, F30, ITSG, Swarm, TOLEOS);
%                                         default none.
%   PoolMonths     1xM   string    -      Whole months (GRACE-FO Level-1B, Swarm, complete SW-All and
%                                         finals2000A), from 2018-06; default none.
%   SwarmMonths    1xS   string    -      Months of Swarm A/B/C densities alone, from 2014-01; default none.
%   DensityMonths  1xD   string    -      Months of TOLEOS densities (CHAMP, GRACE-A, GRACE-FO 1); a
%                                         range or "all" skips the gap 2017-11..2018-04, a single
%                                         month no mission flew in is refused; default none.
%   Items          1xI   string    -      Ids to select (column id); default empty = every id.
%   Today          1x1   datetime  -      The date "ongoing" and "in the future" are judged against;
%                                         default today (UTC). Tests pass a fixed date.
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   T      Rx20  table   -      One row per id: id, group, description, setting (the setup option
%                               that governs it), parameter (the selection that sizes it, "" if fixed),
%                               destination ("toolbox" for the model sources, "pool" for every
%                               data source), unit (what one item is), first, last
%                               (datetime, NaT for undated), ongoing (logical), download_MB, disk_MB
%                               (per item), sizeBasis, licence, licenceText, notice, source (host),
%                               selected (logical), count (items this selection implies; undated
%                               sources count their fixed files), totalDisk_MB (count x disk_MB,
%                               0 when not selected).
%   plan   1x1   struct  -      .PoolDays, .PoolMonths, .SwarmMonths, .DensityMonths (expanded,
%                               1xK string "yyyy-mm-dd" / "yyyy-mm"), .Items (1xI string, the ids
%                               selected), .toleos (1x3 struct: dataset, start, stop as
%                               "yyyy-mm-ddTHH:MM:SS"), .itsgSpans (2x2 datetime, rows GRACE and
%                               GRACE-FO), .swarmFirst, .l1bFirst, .lastComplete, .today (datetime).
%
% ERRORS
%   NonGravForces:fetchCatalog:inputSyntax   An element is not a day, a month, a range or "all".
%   NonGravForces:fetchCatalog:inputRange    A range ends before it starts.
%   NonGravForces:fetchCatalog:inputSpan     A date before the product's first or after 'Today', a
%                                            DensityMonths element that selects no month a TOLEOS
%                                            mission flew in, or one reaching past the last such month.
%   NonGravForces:fetchCatalog:inputItems    An id in Items that the catalogue does not hold.
%
% ACCURACY
%   Sizes are per item and, where sizeBasis says "estimate", good to a factor of two; the measured
%   ones vary by day (an ITSG file 0.34-0.66 MB, an L1B day with the satellites' activity). Counts
%   are exact: they are the number of distinct source URLs setup would request.
%
% REFERENCES
%   THIRD-PARTY-NOTICES.md (licence classes); HAPI 3.0 (startDate/stopDate of the TOLEOS and Swarm
%   headers); HANDOVER Part I (2026-11-06, the catalogue and setup's 'Items').
%
% EXAMPLE
%   [T, plan] = NonGravForces.fetchCatalog('PoolDays', "2020-01-01..2020-01-03");
%   T(T.count > 0, ["id" "count" "totalDisk_MB" "licence"])
%
% See also setupNonGravForces, NonGravForces.dependencyReport.
%
% KEYWORDS fetch catalogue licence data pool installer
% Created 2026-11-06 by Claude Opus 5.5 (Anthropic); every data source to the pool 2026-11-08 by Claude Opus 5.5 (Anthropic); DensityMonths "all" across the TOLEOS gap 2026-11-09 by Claude Opus 5.5 (Anthropic).
arguments
    options.PoolDays (1,:) string = strings(1, 0)
    options.PoolMonths (1,:) string = strings(1, 0)
    options.SwarmMonths (1,:) string = strings(1, 0)
    options.DensityMonths (1,:) string = strings(1, 0)
    options.Items (1,:) string = strings(1, 0)
    options.Today (1,1) datetime = NaT
end
today = options.Today;
if isnat(today)
    today = datetime('now', 'TimeZone', 'UTC');
end
today.TimeZone = '';
today = dateshift(today, 'start', 'day');
lastComplete = dateshift(today, 'start', 'month') - caldays(1);

% THE SPANS (measured; the help says where). Ongoing products end at lastComplete.
itsgSpans = [datetime(2002, 4, 4), datetime(2017, 6, 29); datetime(2018, 6, 1), lastComplete];
swarmFirst = datetime(2014, 1, 1);
l1bFirst = datetime(2018, 6, 1);
toleos = struct('dataset', {"CH_OPER_DNS_ACC_2_", "GR_OPER_DNS1ACC_2_", "GF_OPER_DNS1ACC_2_"}, ...
    'start', {"2000-07-29T00:00:00", "2002-04-04T00:00:00", "2018-05-29T00:00:00"}, ...
    'stop', {"2010-09-04T23:59:50", "2017-10-31T23:59:50", "2026-04-30T23:59:50"});
tStart = datetime([toleos.start], 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');
tStop = datetime([toleos.stop], 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss');

T = localRows(itsgSpans, swarmFirst, l1bFirst, min(tStart), max(tStop), lastComplete);

% Items
items = options.Items;
unknown = setdiff(items, T.id);
if ~isempty(unknown)
    NonGravForces.fail("fetchCatalog", "inputItems", ...
        "unknown id(s) %s -> choose from %s", strjoin(unknown, ", "), strjoin(T.id, ", "));
end
if isempty(items)
    T.selected = true(height(T), 1);
else
    T.selected = ismember(T.id, items);
end

% Selections, expanded
% A reference day has no lower bound (Dst and F30 reach back decades); "all" days are TU Graz's spans.
days = localExpand(options.PoolDays, "day", "PoolDays", itsgSpans, NaT, today);
poolMonths = localExpand(options.PoolMonths, "month", "PoolMonths", [l1bFirst, lastComplete], l1bFirst, today);
swarmMonths = localExpand(options.SwarmMonths, "month", "SwarmMonths", [swarmFirst, lastComplete], swarmFirst, today);
% TOLEOS has a GAP: no mission flew 2017-11..2018-04. "all" is the union of the missions' own spans
% (one row each), never the envelope from the first start to the last stop: that envelope crosses
% the gap, and until 2026-11-09 the check below refused "all" at 2017-11. Each element is expanded
% and judged on its own: it must select at least one flown month and end by the last one; the gap
% months a range crosses are skipped (they hold no file), and the plan lists what remains.
spans = strjoin(string(tStart, 'yyyy-MM') + ".." + string(tStop, 'yyyy-MM'), ", ");
densityMonths = NaT(1, 0);
for e = strtrim(options.DensityMonths)
    em = localExpand(e, "month", "DensityMonths", [tStart(:), tStop(:)], min(tStart), today);
    flown = arrayfun(@(m) any(localFlying(m, tStart, tStop)), em);
    if ~any(flown)
        NonGravForces.fail("fetchCatalog", "inputSpan", ...
            "no TOLEOS mission flew in %s (spans %s) -> pick a month inside one of them", e, spans);
    end
    if any(em(~flown) > dateshift(max(tStop), 'start', 'month'))
        NonGravForces.fail("fetchCatalog", "inputSpan", ...
            "DensityMonths element ""%s"" reaches past the last TOLEOS month %s (spans %s)" + ...
            " -> end it there", e, string(max(tStop), 'yyyy-MM'), spans);
    end
    densityMonths = [densityMonths, em(flown)]; %#ok<AGROW>
end
densityMonths = unique(densityMonths);

% Counts: the distinct source URLs setup requests for this selection, per id
first = unique(dateshift(days, 'start', 'month'));
n = zeros(height(T), 1);
put = @(id) find(T.id == id);
n(put("nrlmsise00")) = 1; n(put("nrlmsis21")) = 1; n(put("hwm14")) = 1; n(put("dtm2020")) = 1;
n(put("eop")) = 1; n(put("spaceweather")) = 1; n(put("igrf")) = 1; n(put("l1b-reference")) = 1;
n(put("ceres")) = 3; n(put("de440")) = 1; n(put("groops-syn1deg")) = 2; n(put("groops-grid")) = 2;
if ~isempty(days)
    n(put("dst")) = numel(unique([first, first - calmonths(1)]));
    n(put("f30")) = 2;
    n(put("itsg")) = 6 * numel(days);          % two satellites x three products per day
end
swarmAll = unique([first(first >= swarmFirst), poolMonths(poolMonths >= swarmFirst), swarmMonths]);
n(put("swarm")) = 3 * numel(swarmAll);
tolMonths = unique([first, densityMonths]);
n(put("toleos")) = sum(arrayfun(@(m) nnz(localFlying(m, tStart, tStop)), tolMonths));
if ~isempty(poolMonths)
    n(put("pool-l1b")) = sum(arrayfun(@(m) day(dateshift(m, 'end', 'month')), poolMonths));
    n(put("pool-spaceweather")) = 1;
    n(put("pool-eop")) = 1;
end
T.count = n;
T.totalDisk_MB = T.count .* T.disk_MB .* T.selected;

plan = struct('PoolDays', string(days, 'yyyy-MM-dd'), 'PoolMonths', string(poolMonths, 'yyyy-MM'), ...
    'SwarmMonths', string(swarmMonths, 'yyyy-MM'), 'DensityMonths', string(densityMonths, 'yyyy-MM'), ...
    'Items', T.id(T.selected)', 'toleos', toleos, 'itsgSpans', itsgSpans, 'swarmFirst', swarmFirst, ...
    'l1bFirst', l1bFirst, 'lastComplete', lastComplete, 'today', today);
end

function T = localRows(itsg, swarmFirst, l1bFirst, tolFirst, tolLast, lastComplete)
%LOCALROWS The catalogue. One row per id; the help says where each span and size was measured.
NaTd = NaT;
oN = "Data files under data/";
r = { ...
% id                 group                         description                                                      setting         parameter                                     dest       unit                first           last          ongoing download disk   sizeBasis                                                        licence      licenceText                                                                                                                                       notice                                                         source
"nrlmsise00",        "Model sources",              "NRLMSISE-00 C source (msis00 density, the default backend)",    "FetchMsis",    "",                                           "toolbox", "git clone",        NaTd,           NaTd,         false,  0.3,     0.2,   "measured on disk 2026-11-06; the clone adds its history",       "unstated",  "The package carries no licence and no redistribution grant, and none is implied by that; it is fetched from its maintainer and never redistributed.", "NRLMSISE-00 — fetched, not redistributed (settled 2026-09-09)", "git.linta.de";
"nrlmsis21",         "Model sources",              "NRLMSIS 2.1 Fortran source (msis21 density, precise preset)",   "FetchMsis",    "",                                           "toolbox", "tarball",          NaTd,           NaTd,         false,  0.3,     0.8,   "measured on disk 2026-11-06",                                   "academic",  "NRL Open Source Academic Research License: research, academic and non-profit use only; commercial use is prohibited.",                      "NRLMSIS 2.1 — fetched, not redistributed (settled 2026-09-15)", "map.nrl.navy.mil";
"hwm14",             "Model sources",              "HWM14 Fortran source and data (drag.wind hwm14, opt-in)",          "FetchMsis",    "",                                           "toolbox", "7 files",          NaTd,           NaTd,         false,  0.3,     0.3,   "measured on disk 2026-11-06",                                   "unstated",  "NRL's package carries no licence or redistribution statement; two data files come from the MIT-licensed pyHWM14 repository.",                "HWM14 — fetched, not redistributed (settled 2026-09-15)",     "map.nrl.navy.mil";
"dtm2020",           "Model sources",              "DTM2020 Fortran source and coefficients (opt-in density)",      "FetchMsis",    "",                                           "toolbox", "3 files",          NaTd,           NaTd,         false,  1,       1,     "estimate",                                                      "academic",  "CNES: academic, non-commercial use only; no modification, no redistribution; name CNES's DTM in any publication based on it.",               "DTM2020 — fetched, not redistributed (2026-10-27)",           "raw.githubusercontent.com";
"eop",               "Reference data",             "IERS finals2000A.all Earth orientation",                        "FetchData",    "",                                           "pool",    "file",             NaTd,           NaTd,         false,  3.8,     3.8,   "measured 2026-11-06",                                           "unreviewed", "IERS product, under the provider's terms; this repository has not reviewed them.",                                                         "",                                                            "datacenter.iers.org";
"spaceweather",      "Reference data",             "CelesTrak SW-All.csv space weather",                            "FetchData",    "",                                           "pool",    "file",             NaTd,           NaTd,         false,  2.9,     2.9,   "measured 2026-11-06",                                           "unreviewed", "CelesTrak product, under the provider's terms; this repository has not reviewed them.",                                                     "",                                                            "celestrak.org";
"igrf",              "Reference data",             "IGRF-14 coefficients (IAGA V-MOD)",                             "FetchData",    "",                                           "pool",    "file",             NaTd,           NaTd,         false,  0.04,    0.04,  "measured on disk 2026-11-06",                                   "unreviewed", "IAGA V-MOD product, under the provider's terms; this repository has not reviewed them.",                                                     "",                                                            "www.ngdc.noaa.gov";
"l1b-reference",     "Reference data",             "GRACE-FO Level-1B GNV1B/SCA1B/MAG1B of 2020-01-01 (C)",         "FetchData",    "",                                           "pool",    "day archive",      NaTd,           NaTd,         false,  147,     40,    "download as stated by setup; disk estimate",                    "unreviewed", "NASA/JPL GRACE-FO Level-1B served by GFZ's ISDC, under the provider's terms; this repository has not reviewed them.",                       "",                                                            "isdc-data.gfz.de";
"ceres",             "Reference data",             "CERES EBAF 4.2.1 monthly means rsut/rlut/rsdt 2000-03..2025-09", "FetchData",   "",                                           "pool",    "file",             NaTd,           NaTd,         false,  21.7,    21.7,  "65 MB for the three, as stated by setup",                       "open",      "NASA LaRC, CC BY 4.0: reuse allowed, commercial included, with attribution (data/ceres/README.md).",                                        oN,                                                            "esgf-node.ornl.gov";
"de440",             "Ephemeris",                  "JPL DE440 planetary ephemeris (de440s.bsp)",                    "FetchKernels", "",                                           "pool",    "file",             NaTd,           NaTd,         false,  32.7,    32.7,  "measured: 32 726 016 bytes",                                    "unreviewed", "NASA/JPL NAIF kernel, under the provider's terms; this repository has not reviewed them.",                                                  "",                                                            "naif.jpl.nasa.gov";
"groops-syn1deg",    "Earth-radiation grids",      "GROOPS CERES SYN1deg albedo and longwave climatologies",       "FetchGroops",  "",                                           "pool",    "file",             NaTd,           NaTd,         false,  1250.9,  1250.9, "measured in the pool 2026-11-06",                              "unreviewed", "Published openly in TU Graz's GROOPS data pool, derived from NASA CERES SYN1deg; the terms are not recorded here.",                          "",                                                            "ftp.tugraz.at";
"groops-grid",       "Earth-radiation grids",      "GROOPS CERES reflectivity and emissivity 2.5 deg grids",        "FetchGroops",  "",                                           "pool",    "file",             NaTd,           NaTd,         false,  1.1,     1.1,   "measured in the pool 2026-11-06",                               "unreviewed", "Published openly in TU Graz's GROOPS data pool, derived from NASA CERES; the terms are not recorded here.",                                  "",                                                            "ftp.tugraz.at";
"dst",               "Reference-day observations", "Kyoto WDC final hourly Dst (each day's month and the one before)", "FetchPool",  "PoolDays",                                   "pool",    "month",            NaTd,           NaTd,         false,  0.03,    0.03,  "estimate",                                                      "unreviewed", "World Data Center for Geomagnetism, Kyoto, under the provider's terms; this repository has not reviewed them.",                              "",                                                            "wdc.kugi.kyoto-u.ac.jp";
"f30",               "Reference-day observations", "CLS 30 cm radio flux F30 (adjusted and absolute files)",        "FetchPool",    "PoolDays",                                   "pool",    "file",             NaTd,           NaTd,         false,  1,       1,     "estimate",                                                      "unreviewed", "CLS product, under the provider's terms; this repository has not reviewed them.",                                                          "",                                                            "ftpsedr.cls.fr";
"itsg",              "Reference-day observations", "TU Graz attitude, non-conservative forces and orbits (both satellites)", "FetchPool", "PoolDays",                           "pool",    "file",             itsg(1, 1),     itsg(2, 2),   true,   0.48,    0.48,  "measured 2026-11-06: 0.34-0.66 MB per file",                    "unreviewed", "TU Graz ITSG operational products, under the provider's terms; this repository has not reviewed them.",                                      "",                                                            "ftp.tugraz.at";
"swarm",             "Densities",                  "ESA Swarm A/B/C POD densities DNSxPOD (one file per satellite and month)", "FetchPool", "PoolDays, PoolMonths, SwarmMonths",  "pool",    "file",             swarmFirst,     lastComplete, true,   13.3,    13.3,  "about 40 MB per month for three, as stated by setup",           "open",      "ESA Swarm data, free and open under ESA's terms (vires.services/data_terms); credit: Data provided by ESA, Swarm mission.",               "ESA Swarm thermosphere densities (2026-10-27)",              "vires.services";
"toleos",            "Densities",                  "TU Delft accelerometer densities, ESA TOLEOS (CHAMP, GRACE-A, GRACE-FO 1)", "FetchPool", "PoolDays, DensityMonths",          "pool",    "file",             tolFirst,       tolLast,      false,  39.9,    39.9,  "measured 2026-11-06: a GRACE-FO month",                         "open",      "ESA TOLEOS version 02, CC BY 4.0: reuse allowed, commercial included, with attribution to ESA and TU Delft.",                                 "ESA TOLEOS accelerometer densities (2026-10-28)",            "vires.services";
"pool-l1b",          "Whole months",               "GRACE-FO Level-1B ACT1B/GNV1B/SCA1B/THR1B/IMU1B of C and D, every day", "FetchPool", "PoolMonths",                             "pool",    "day archive",      l1bFirst,       lastComplete, true,   150,     299,   "download as stated by setup; disk measured 2026-11-06 (34 days)", "unreviewed", "NASA/JPL GRACE-FO Level-1B served by GFZ's ISDC, under the provider's terms; this repository has not reviewed them.",                     "",                                                            "isdc-data.gfz.de";
"pool-spaceweather", "Whole months",               "CelesTrak SW-All.csv, complete, into the pool",                 "FetchPool",    "PoolMonths",                                 "pool",    "file",             NaTd,           NaTd,         false,  2.9,     2.9,   "measured 2026-11-06",                                           "unreviewed", "CelesTrak product, under the provider's terms; this repository has not reviewed them.",                                                     "",                                                            "celestrak.org";
"pool-eop",          "Whole months",               "IERS finals2000A.all, complete, into the pool",                 "FetchPool",    "PoolMonths",                                 "pool",    "file",             NaTd,           NaTd,         false,  3.8,     3.8,   "measured 2026-11-06",                                           "unreviewed", "IERS product, under the provider's terms; this repository has not reviewed them.",                                                         "",                                                            "datacenter.iers.org"};
names = ["id" "group" "description" "setting" "parameter" "destination" "unit" "first" "last" "ongoing" ...
    "download_MB" "disk_MB" "sizeBasis" "licence" "licenceText" "notice" "source"];
T = cell2table(r, 'VariableNames', names);
for v = ["id" "group" "description" "setting" "parameter" "destination" "unit" "sizeBasis" "licence" "licenceText" "notice" "source"]
    T.(v) = string(T.(v));
end
end

function d = localExpand(list, kind, name, spans, lower, today)
%LOCALEXPAND Days or months from singles, "a..b" ranges and "all"; sorted and unique (datetime row).
%   spans (Sx2) is what "all" means; lower (NaT = none) and today bound an explicit element.
d = NaT(1, 0);
if kind == "day"
    pat = digitsPattern(4) + "-" + digitsPattern(2) + "-" + digitsPattern(2);
    fmt = 'yyyy-MM-dd';
else
    pat = digitsPattern(4) + "-" + digitsPattern(2);
    fmt = 'yyyy-MM';
end
for e = strtrim(list)
    if e == "all"
        for s = 1:size(spans, 1)
            a = spans(s, 1); b = spans(s, 2);
            if kind == "month", a = dateshift(a, 'start', 'month'); b = dateshift(b, 'start', 'month'); end
            d = [d, localStep(a, b, kind)]; %#ok<AGROW>
        end
        continue
    end
    parts = split(e, "..");
    if ~(numel(parts) <= 2 && all(matches(parts, pat)))
        NonGravForces.fail("fetchCatalog", "inputSyntax", ...
            "%s element ""%s"" is not ""%s"", ""%s..%s"" or ""all"" -> e.g. ""%s""", ...
            name, e, fmt, fmt, fmt, string(today, fmt));
    end
    % A date that does not round-trip is not a calendar date ("2021-02-29" would roll over to March).
    try
        a = datetime(parts(1), 'InputFormat', fmt);
        b = datetime(parts(end), 'InputFormat', fmt);
        ok = ~any(isnat([a b])) && string(a, fmt) == parts(1) && string(b, fmt) == parts(end);
    catch
        ok = false;
    end
    if ~ok
        NonGravForces.fail("fetchCatalog", "inputSyntax", "%s element ""%s"" is not a calendar date -> correct it", name, e);
    end
    if b < a
        NonGravForces.fail("fetchCatalog", "inputRange", "%s range ""%s"" ends before it starts -> swap its ends", name, e);
    end
    first = lower;
    if kind == "month" && ~isnat(first), first = dateshift(first, 'start', 'month'); end
    if (~isnat(first) && a < first) || b > today
        if isnat(first), from = "any date"; else, from = string(first, fmt); end
        NonGravForces.fail("fetchCatalog", "inputSpan", ...
            "%s element ""%s"" reaches outside %s..%s -> stay inside the published span", ...
            name, e, from, string(today, fmt));
    end
    d = [d, localStep(a, b, kind)]; %#ok<AGROW>
end
d = unique(d);
end

function d = localStep(a, b, kind)
if kind == "day", d = a:caldays(1):b; else, d = a:calmonths(1):b; end
end

function tf = localFlying(m, tStart, tStop)
%LOCALFLYING Which TOLEOS datasets cover part of month m; setup's localToleosTargets clips the same way.
t0 = max(m, tStart);
t1 = min(m + calmonths(1), tStop);
tf = t1 > t0;
end
