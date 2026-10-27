function info = setupNonGravForces(options)
%SETUPNONGRAVFORCES Put the NonGravForces toolbox on the path, fetch what it does not ship, report the dependencies.
%
%   setupNonGravForces
%   info = setupNonGravForces('Verbose', false, 'FetchData', "never")
%
%   Adds the toolbox root, tests/, tools/ and demos/ to the MATLAB path (no savepath),
%   registers the documentation for the Help browser, and reports which of the
%   ecosystem dependencies are reachable. It never adds a dependency to the path
%   itself: each dependency ships its own setup function and owns its own path.
%
%   IT IS THE ONE PLACE IN THE TOOLBOX THAT TOUCHES THE NETWORK. What it fetches,
%   and why it is fetched rather than shipped:
%
%     model sources ('FetchMsis')
%       NRLMSISE-00 C       git clone at a pinned commit into mex/nrlmsise00 (no
%                           redistribution grant, O-34)
%       NRLMSIS 2.1 Fortran NRL's tarball, pinned by SHA-256, unpacked into
%                           mex/nrlmsis21 (Open Source Academic Research License:
%                           research and non-profit use only; the statement is printed
%                           when it is fetched and recorded in THIRD-PARTY-NOTICES.md)
%       HWM14 Fortran       NRL's release files, each pinned by SHA-256, into mex/hwm14; the
%                           two DWM/QD data files NRL's server refuses (HTTP 403, 2026-09-15)
%                           from the pyHWM14 repository at a pinned commit, where the model's
%                           own check output, reproduced exactly, is what vouches for them
%       DTM2020 Fortran     CNES's operational DTM2020 (F10.7, Kp) as SWAMI's MCM repository
%                           ships it, at a pinned commit, each file pinned by SHA-256, into
%                           mex/dtm2020: the subroutine file, the coefficient file and the
%                           licence (CNES: academic, non-commercial use only; no modification,
%                           no redistribution -- printed when fetched, THIRD-PARTY-NOTICES.md)
%     data ('FetchData'), the DE440 ephemeris ('FetchKernels'), the GROOPS Earth-radiation
%     grids ('FetchGroops', into DATAPOOL)
%       finals2000A.all     IERS Earth orientation, into data/eop
%       SW-All.csv          CelesTrak space weather, into data/spaceweather
%       igrf14coeffs.txt    IAGA V-MOD IGRF-14 coefficients, into data/igrf
%       GRACE-FO Level-1B   GNV1B, SCA1B and MAG1B of 2020-01-01 (GRACE-C), full day, from the
%                           archive of GFZ's ISDC, into data/l1b -- NOT data/samples, whose
%                           committed ten-minute excerpts carry the same file names and are
%                           the tests' oracles
%       CERES EBAF 4.2.1    rsut, rlut and rsdt monthly means 2000-03..2025-09 (obs4MIPs,
%                           65 MB) from ESGF, each kept only if its SHA-256 is the one ESGF
%                           publishes, into data/ceres beside the committed two-month excerpt
%     observations for the reference days ('FetchPool', into DATAPOOL; PR #107)
%       Dst                 Kyoto WDC FINAL hourly Dst, the month of each 'PoolDays' day and the
%                           month before it, into <DataPool>/SpaceWeather/Dst
%       F30                 CLS daily 30 cm radio flux (adjusted to 1 AU and absolute; gaps
%                           interpolated and flagged), by anonymous FTP, into
%                           <DataPool>/SpaceWeather/F30
%       TU Graz products    attitude, nonConservativeForces and reducedDynamicOrbit of GRACE-A/B
%                           (days before 2018) or GRACE-C/D (from 2018) per 'PoolDays' day, as
%                           published (.txt.gz), into <DataPool>/ITSG/operational/<sat>/<product>/<yyyy>
%       Swarm densities     ESA's Level-2 SW_OPER_DNSxPOD_2_ (thermosphere density from precise orbit
%                           determination alone, 30 s) of Swarm A, B and C for the month of each
%                           'PoolDays' day from 2014 on, through the VirES for Swarm HAPI server as
%                           CSV with its JSON header, into <DataPool>/Swarm/DNSxPOD/<yyyy>
%                           The committed excerpts in data/spaceweather and data/samples/itsg are
%                           cut from these by tools/cutFixtures.py. Classic GRACE Level-1B is not
%                           here yet: GFZ's ISDC, its only anonymous source, is down (2026-10-26)
%                           and no archive name has been read from it.
%     whole months ('PoolMonths', into DATAPOOL, with FetchPool's setting; O-19)
%       GRACE-FO Level-1B   ACT1B, GNV1B, SCA1B and THR1B of GRACE-C and GRACE-D for every day of
%                           each month, members only, from ISDC's daily archives into
%                           <DataPool>/GRACEFO/L1B/RL04/<yyyy>/<yyyy-mm-dd> (tools/dragScaleMonth.m)
%       Swarm densities     as above, for each month
%       SW-All, finals2000A CelesTrak's complete space-weather file and IERS finals2000A.all into
%                           <DataPool>/SpaceWeather and <DataPool>/EOP; refetched when the file in the
%                           pool ends before the last month (a stale copy is replaced, not trusted)
%
%   PRESENCE IS DECIDED BY THE EXACT FILE, AND A DOWNLOAD BY ITS CONTENT. Until
%   2026-09-15 a glob decided presence, so the committed seven-month excerpts
%   (finals2000A_2019-09_2020-03.all, ...) counted as the full products and the
%   full products were never fetched; and whatever a server returned was kept, an
%   HTML error page included (HANDOVER LP-C, finding F15; CLAUDE.md: "a status code
%   is not a file"). Every fetched file is now checked by its header or digest
%   before it replaces anything, and a failed check is REPORTED, never raised.
%   The Level-1B members are read out of the archive by this function's own tar
%   reader: MATLAB's untar takes no member list (its three-argument call, used
%   until 2026-09-15, raised MATLAB:narginchk:tooManyInputs on every attempt, and
%   the failure was reported as a download failure).
%
%   A build step never fetches: NonGravForces.buildMex refuses when a source is
%   missing and names this function, because a build that reaches the network
%   without being asked is a build nobody can audit. Every network action prints
%   one line whatever 'Verbose' says.
%
% INPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   (none)
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type     Units  Description
%   ----       ----  ----     -----  -----------
%   Verbose    1x1   logical  -      Print the dependency and fetch table (default true).
%   FetchMsis  1x1   string   -      "ifMissing" (default) fetches the model sources only when absent;
%                                    "never" leaves the network alone; "force" fetches again.
%   FetchData  1x1   string   -      The same three settings for the data files.
%   FetchGroops 1x1 string    -      The same three settings for the GROOPS Earth-radiation data,
%                                    DEFAULT "never": the two CERES SYN1deg climatologies are
%                                    1.25 GB EACH, so they are never pulled unasked. They are not
%                                    committed either (O-41); they live in DATAPOOL.
%   DataPool   1x1   string   -      Root of the GROOPS-style data pool. Default: the DATAPOOL
%                                    environment variable. It is NOT guessed: a 2.5 GB download
%                                    into a directory nobody named is not a service (R5). The
%                                    GROOPS files go into <DataPool>/EarthRadiationPressure,
%                                    FetchPool's into <DataPool>/SpaceWeather, <DataPool>/ITSG and <DataPool>/Swarm.
%   FetchPool  1x1   string   -      The same three settings for the observations of the reference
%                                    days (Dst, F30, TU Graz products, Swarm densities), into DataPool (default
%                                    "ifMissing"; nothing happens without a DataPool). ciGate passes
%                                    "never": the gate reads committed excerpts, not the network.
%   PoolDays   1xN   string   -      Reference days for FetchPool, "yyyy-mm-dd" (default the fixture
%                                    days "2008-01-01" "2020-01-01" "2020-01-02" "2020-01-03").
%   PoolMonths 1xM   string   -      Whole months for FetchPool, "yyyy-mm" from 2018-06 (default none):
%                                    GRACE-C/D Level-1B, Swarm densities and the complete space-weather
%                                    and EOP files a month run (tools/dragScaleMonth.m) needs. A day
%                                    archive is about 150 MB, a month about 4.5 GB.
%   FetchKernels 1x1 string   -      The same three settings for the DE440 ephemeris. It is SEPARATE
%                                    from FetchData on purpose: ciGate passes FetchData "never" so
%                                    that the gate does not test the network for OBSERVATIONS, which
%                                    are committed as excerpts. A planetary ephemeris is model data,
%                                    like the MSIS and HWM14 sources, and the gate needs it or the
%                                    precise preset's own sun.source stays untested (HANDOVER O-7).
%
% OUTPUT ARGUMENTS
%   Name      Size  Type     Units  Description
%   ----      ----  ----     -----  -----------
%   info      1x1   struct   -      .root (char), .version (char),
%                                   .dependencies (struct array: name, required, present, version, setup),
%                                   .mex (1xT NonGravForces.buildMex('Check') report, one element per target),
%                                   .msis (.dir, .action "present"|"fetched"|"skipped"|"failed", .url, .commit, .message),
%                                   .sources (struct array: name, action, dir, source, message),
%                                   .data (.items struct array: name, action, file, source; .mode;
%                                   the Level-1B item's file lists the members, one per line)
%
% ERRORS
%   None raised for absent dependencies or failed fetches: absence is REPORTED (rule R4/R5).
%
% EXAMPLE
%   info = setupNonGravForces('FetchData', "never");
%   [info.sources.name; info.sources.action]
%
% See also NonGravForces.dependencyReport, NonGravForces.buildMex, runNgfTests.
%
% KEYWORDS setup path dependencies fetch
% Created 2026-08-27 by Claude Fable 5 (Anthropic) for Matthias Weigelt; content-checked fetches and NRLMSIS 2.1 2026-09-15 by Claude Opus 5 (Anthropic); FetchPool 2026-10-26 by Claude Opus 5.5 (Anthropic).
arguments
    options.Verbose (1,1) logical = true
    options.FetchMsis (1,1) string {mustBeMember(options.FetchMsis, ["ifMissing","never","force"])} = "ifMissing"
    options.FetchData (1,1) string {mustBeMember(options.FetchData, ["ifMissing","never","force"])} = "ifMissing"
    options.FetchKernels (1,1) string {mustBeMember(options.FetchKernels, ["ifMissing","never","force"])} = "ifMissing"
    options.FetchGroops (1,1) string {mustBeMember(options.FetchGroops, ["ifMissing","never","force"])} = "never"
    options.DataPool (1,1) string = string(getenv('DATAPOOL'))
    options.FetchPool (1,1) string {mustBeMember(options.FetchPool, ["ifMissing","never","force"])} = "ifMissing"
    options.PoolDays (1,:) string {mustBeNonempty} = ["2008-01-01" "2020-01-01" "2020-01-02" "2020-01-03"]
    options.PoolMonths (1,:) string {localMustBeMonths} = strings(1, 0)
end
% THE PINS AND URLS LIVE HERE AND NOWHERE ELSE (rule R1). Moving a pin is an act with a run behind it.
% NRLMSISE-00: bc9a2fe is "gtd7: fix output at 32.5 km" (Argon densities and a memory corruption fixed
% since the 2007 release this toolbox used until 2026-09-09).
MSIS_URL = "https://git.linta.de/~brodo/nrlmsise-00.git";
MSIS_COMMIT = "bc9a2feba4344e74201281e563332688a4d09cc3";
% NRLMSIS 2.1: NRL's release tarball of 2022-04-04, digest measured 2026-09-15 when it was first fetched.
MSIS21_URL = "https://map.nrl.navy.mil/map/pub/nrl/NRLMSIS/NRLMSIS2.1/nrlmsis2.1.tar.gz";
MSIS21_SHA256 = "41e47b29f795d36a5cc252b2858aa2a384c4a7323ace3d48d3ea2f2b37a1a6a8";
% HWM14: NRL's supplement to Drob et al. (2015), and the two data files NRL's server answers 403 for,
% from rilma/pyHWM14 at 115fd958 (its hwm123114.bin is byte-identical to NRL's). Digests measured 2026-09-15,
% when checkhwm14 built from exactly these files reproduced NRL's Check/gfortran.txt line for line.
HWM14_NRL = "https://map.nrl.navy.mil/map/pub/nrl/HWM/HWM14/HWM14_ess224-sup-0002-supinfo/";
HWM14_PY = "https://raw.githubusercontent.com/rilma/pyHWM14/115fd958f4d196d6eb89bf6d997b5930355ada98/data/";
HWM14_FILES = struct( ...
    'name', {"hwm14.f90", "hwm123114.bin", "checkhwm14.f90", "gfortran.txt", "README.txt", "dwm07b104i.dat", "gd2qd.dat"}, ...
    'url',  {HWM14_NRL + "hwm14.f90", HWM14_NRL + "hwm123114.bin", HWM14_NRL + "checkhwm14.f90", HWM14_NRL + "Check/gfortran.txt", ...
             HWM14_NRL + "README.txt", HWM14_PY + "dwm07b104i.dat", HWM14_PY + "gd2qd.dat"}, ...
    'sha256', {"6bb4e031917b44c93201f0289ada8bf6f6d81eb3f292c89151c6937e2e7bf76e", "6e445f8337c7efc815b7ff7f9f967d8a9a16b469d93df9c46e54930553beb906", ...
               "6a58eaafb065ffa9d34cfba409f865531b8cbe7c1a0ee274ca6c21fcbeb26091", "2b1d4f4f103be3531393c48bf32d034548babfb6c080549884cad9fd3a2c8652", ...
               "14b6e5e48d346ff034e2f6bfd7f26dc550945c9b90ad9661c7710e721ac89b15", "f2b8eff002d55b0f6d49d202c73b7a9cb6685f2c6bb57f1291a2adc7aa07cf4f", ...
               "6bb1f2384e30b409240ee92c32726ed2a3d73eb550c05e0969d1178032319804"});
% DTM2020 (operational: F10.7 and Kp; Bruinsma & Boniface 2021), CNES, as distributed by the SWAMI project.
% The libswamif subroutine file is used because its lecdtm takes the unit number the documentation asks
% the caller to define; the src/dtm2020 copy opens a fixed file name and reads an undefined unit, and the
% licence forbids editing it. Above 150 km the two dtm3 are the same code (the MCM file adds a blend
% below 150 km for its whole-atmosphere model). Digests measured 2026-10-27 from git blobs at the pin.
DTM2020_BASE = "https://raw.githubusercontent.com/swami-h2020-eu/mcm/a488a7c9d030bfbe86e88ab3d28a7ec5589b92e0/";
DTM2020_FILES = struct( ...
    'name', {"dtm2020_F107_Kp-subr_MCM.f90", "DTM_2020_F107_Kp.dat", "LICENSE_MCM_DTM.txt"}, ...
    'url',  {DTM2020_BASE + "src/libswamif/dtm2020_F107_Kp-subr_MCM.f90", DTM2020_BASE + "data/DTM_2020_F107_Kp.dat", DTM2020_BASE + "LICENSE"}, ...
    'sha256', {"1b109af61c85b90177d028e89fa0a6ba1c7208fb6e0f50fb1a793062a1e2f208", "c08f7ca61a3521a736a857850bd6cacd169f9302c356dc8523608b1a1d941588", ...
               "b45cc69494a849efd7232d3e2de85480c8ac6e389c974244a9b882aeab4e5085"});
EOP_URL = "https://datacenter.iers.org/data/9/finals2000A.all";
SW_URL = "https://celestrak.org/SpaceData/SW-All.csv";
% IGRF-14, IAGA V-MOD. Verified 2026-09-09: byte-identical to the SHA-256 data/igrf/README.md records.
IGRF_URL = "https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt";
% GRACE-FO Level-1B is anonymously available at GFZ's ISDC (owner, 2026-09-09): one archive per day.
L1B_BASE = "https://isdc-data.gfz.de/grace-fo/Level-1B/JPL/INSTRUMENT/RL04";
L1B_DAY = "2020-01-01";
L1B_MEMBERS = ["GNV1B_2020-01-01_C_04.txt", "SCA1B_2020-01-01_C_04.txt", "MAG1B_2020-01-01_C_04.txt"];
% What tools/dragScaleMonth.m reads of a month day, for both satellites (FetchPool 'PoolMonths').
L1B_MONTH_PRODUCTS = ["ACT1B" "GNV1B" "SCA1B" "THR1B"];
% CERES EBAF 4.2.1 on ESGF (obs4MIPs, dataset version v20260220), anonymous; the SHA-256 are ESGF's own
% (esgf-1-5-bridge search, 2026-09-15), and each download was checked against them the same day.
CERES_BASE = "https://esgf-node.ornl.gov/thredds/fileServer/user_pub_work/obs4MIPs/NASA-LaRC/CERES-EBAF-4-2-1/mon";
CERES_FILES = ["rsut" "rlut" "rsdt"] + "_mon_CERES-EBAF-4-2-1_NASA-GISS_gn_200003-202509.nc";
CERES_SHA256 = ["6869b2207d86ec6581a6518f63b65f41fc209dda4ea5e64afbf6f3b6d47ea592", ...
                "1e611d51416d512a945ee4a61983f5c88a991f08a6aee6d6d7662ea76bd9ab5a", ...
                "3e4cd4ad8f1c283369cd2dd385ac3f655c4e10c9ed2a6b954b829097a61125e8"];
% DE440: the planetary ephemeris the precise preset names. THE URL IS NOT HERE -- 
% referenceSystem.downloadKernels calls itself "the SINGLE AUTHORITY for where these files come
% from" and it is, so this toolbox asks that catalogue to fetch and pins only WHAT ARRIVES.
% de440s.bsp is the short issue (1849-2150), ample for anything GRACE-era, 32 726 016 bytes.
% THE DIGEST IS MEASURED HERE, NOT QUOTED FROM THE FETCHER. The first value written into this
% line was the one downloadKernels REPORTS in its manifest, and it describes no file this project
% can find: MATLAB's own java digest, Python's hashlib and the CI runner all say c1c7feea... for
% the bytes on disk. Taking a number from the tool whose work you are checking is not a check
% (upstream finding, HANDOVER O-40).
% GROOPS publishes its data pool openly at TU Graz (owner, 2026-09-17), which is what closes the
% redistribution question these files raised: they are FETCHED, not vendored. The two SYN1deg
% climatologies are 1.25 GB each and live in DATAPOOL, never in this repository.
GROOPS_BASE = "https://ftp.tugraz.at/pub/ITSG/groops/data";
GROOPS_FILES = struct( ...
    'name', {"CERES_SYN1deg_albedo_climatology.dat", "CERES_SYN1deg_longwaveFlux_climatology.dat", ...
             "earth_ceres_reflectivity_grid2.5.dat", "earth_ceres_emissivity_grid2.5.dat"}, ...
    'folder', {"earthRadiationPressure", "earthRadiationPressure", "albedo", "albedo"});
% Observations for the reference days (FetchPool). The three sources were found by the fetch-data
% workflow's probes of 2026-10-26 (PR #107, runs 40358650063 to 40375328241). NO SHA-256 PINS: the
% days are a parameter, and CLS rewrites its files daily; each reply is checked by its first line.
DST_BASE = "https://wdc.kugi.kyoto-u.ac.jp/dst_final";      % /<yyyymm>/dst<yymm>.for.request
F30_BASE = "ftp://ftpsedr.cls.fr/pub/previsol/solarflux/observation";
F30_FILES = ["radio_flux_adjusted_observation.txt" "radio_flux_absolute_observation.txt"];
ITSG_BASE = "https://ftp.tugraz.at/pub/ITSG/satelliteOrbitProducts/operational";
ITSG_PRODUCTS = ["attitude" "nonConservativeForces" "reducedDynamicOrbit"];
% Swarm's POD-only densities (ESA Level 2, free and open under ESA's data policy): an observed
% density independent of GRACE-FO's accelerometer, for the drag scale of O-19 (2026-10-27).
% HAPI 3 request syntax; a month per satellite per request, the header included so the file
% describes itself and the content check has a first line to read.
SWARM_HAPI = "https://vires.services/hapi/data";
SWARM_SATS = ["A" "B" "C"];
DE440_FILE = "de440s.bsp";
DE440_SHA256 = "c1c7feeab882263fc493a9d5a5b2ddd71b54826cdf65d8d17a76126b260a49f2";
root = fileparts(mfilename('fullpath'));
addpath(root, fullfile(root,'tests'), fullfile(root,'tools'), fullfile(root,'demos'));
% Help-browser registration: info.xml at the root is picked up by doc.
try
    builddocsearchdb(fullfile(root,'doc','html'));
catch
    % Search database is optional; the pages are readable without it.
end
info.root    = root;
info.version = NonGravForces.toolboxVersion('string');
info.dependencies = NonGravForces.dependencyReport('Verbose', options.Verbose);
info.msis = localFetchMsis(root, MSIS_URL, MSIS_COMMIT, options.FetchMsis);
msis21 = localFetchArchive("nrlmsis21", fullfile(root, 'mex', 'nrlmsis21'), "msis_calc.F90", MSIS21_URL, MSIS21_SHA256, ...
    options.FetchMsis, "This software incorporates the MSIS(R) empirical atmospheric model software designed and provided by NRL. " + ...
    "Use is governed by the Open Source Academic Research License Agreement contained in the file nrlmsis2.1_license.txt (research, academic and non-profit use only).");
hwm14 = localFetchFiles("hwm14", fullfile(root, 'mex', 'hwm14'), HWM14_FILES, options.FetchMsis);
dtm2020 = localFetchFiles("dtm2020", fullfile(root, 'mex', 'dtm2020'), DTM2020_FILES, options.FetchMsis);
if dtm2020.action == "fetched"
    dtm2020.message = "The DTM software is property of CNES, licensed for academic, non-commercial purposes only; it may not be modified " + ...
        "or redistributed (mex/dtm2020/LICENSE_MCM_DTM.txt). Mention CNES's DTM in any publication of results based on it.";
    fprintf('setupNonGravForces: %s\n', dtm2020.message);
end
info.sources = [struct('name', "nrlmsise00", 'action', info.msis.action, 'dir', info.msis.dir, 'source', MSIS_URL, 'message', info.msis.message), msis21, hwm14, dtm2020];
info.data = localFetchData(root, EOP_URL, SW_URL, IGRF_URL, L1B_BASE, L1B_DAY, L1B_MEMBERS, options.FetchData);
info.data.items = [info.data.items, localFetchCeres(root, CERES_BASE, CERES_FILES, CERES_SHA256, options.FetchData)];
info.data.items = [info.data.items, localFetchKernel(root, DE440_FILE, DE440_SHA256, options.FetchKernels)];
info.dataPool = options.DataPool;
info.data.items = [info.data.items, localFetchGroops(options.DataPool, GROOPS_BASE, GROOPS_FILES, options.FetchGroops)];
info.data.items = [info.data.items, localFetchPool(options.DataPool, options.PoolDays, DST_BASE, F30_BASE, F30_FILES, ...
    ITSG_BASE, ITSG_PRODUCTS, SWARM_HAPI, SWARM_SATS, options.FetchPool)];
info.data.items = [info.data.items, localFetchPoolMonths(options.DataPool, options.PoolMonths, L1B_BASE, L1B_MONTH_PRODUCTS, ...
    SWARM_HAPI, SWARM_SATS, EOP_URL, SW_URL, options.FetchPool)];
info.mex = NonGravForces.buildMex('Check', true);      % 1xT: one element per target
if options.Verbose
    for m = 1:numel(info.mex)
        if info.mex(m).upToDate, state = 'built, current';
        elseif info.mex(m).present, state = 'built, STALE (sources changed): run NonGravForces.buildMex';
        else, state = 'NOT BUILT: run NonGravForces.buildMex (the caller of this gateway refuses until then)';
        end
        [~, gw] = fileparts(info.mex(m).file);
        fprintf('MEX  %-14s %s\n', gw, state);
    end
    for s = 1:numel(info.sources)
        fprintf('SRC  %-14s %s (%s)\n', info.sources(s).name, info.sources(s).action, info.sources(s).dir);
    end
    for d = 1:numel(info.data.items)
        fprintf('DATA %-14s %s\n', info.data.items(d).name, info.data.items(d).action);
    end
end
end

function d = localFetchData(root, eopUrl, swUrl, igrfUrl, l1bBase, l1bDay, l1bMembers, mode)
%LOCALFETCHDATA The data files that can be had without an account, each checked by its content.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
% name, folder, exact file, url, content check (a pattern the file's head must match)
targets = {"eop", fullfile(root, 'data', 'eop'), "finals2000A.all", eopUrl, '^\s*\d{1,2}\s*\d{1,2}\s*\d{1,2}\s+\d{5}\.\d{2}\s+[IP]'; ...
           "spaceweather", fullfile(root, 'data', 'spaceweather'), "SW-All.csv", swUrl, '^DATE,BSRN,ND,KP1'; ...
           "igrf", fullfile(root, 'data', 'igrf'), "igrf14coeffs.txt", igrfUrl, 'IGRF'};
for t = 1:size(targets, 1)
    folder = targets{t, 2};
    dst = fullfile(folder, char(targets{t, 3}));
    present = isfile(dst);
    rec = struct('name', targets{t, 1}, 'action', "present", 'file', string(dst), 'source', targets{t, 4});
    if mode == "never" || (present && mode ~= "force")
        if ~present, rec.action = "missing (FetchData never)"; rec.file = ""; end
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if ~isfolder(folder), mkdir(folder); end
    tmp = [tempname '.download'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching %s from %s\n', targets{t, 1}, targets{t, 4});
    ok = false;
    try
        localDownload(tmp, targets{t, 4});  % CelesTrak answered 503 once and 200 on the retry
        ok = true;
    catch err
        rec.action = "failed: " + string(err.message);
    end
    if ok && ~localHeadMatches(tmp, targets{t, 5})
        ok = false;
        rec.action = "failed: the server's reply is not the expected product (its first lines do not match) -> it was discarded";
    end
    if ok
        movefile(tmp, dst, 'f');
        rec.action = "fetched";
    else
        rec.file = "";
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, targets{t, 4});
    end
    items(end+1) = rec; %#ok<AGROW>
    clear cleanup
end
% GRACE-FO Level-1B: one archive per day, 147 MB, fetched only when a member is absent, members only.
l1b = fullfile(root, 'data', 'l1b');
url = l1bBase + "/" + extractBefore(l1bDay, 5) + "/gracefo_1B_" + l1bDay + "_RL04.ascii.noLRI.tgz";
missing = l1bMembers(~isfile(fullfile(l1b, l1bMembers)));
rec = struct('name', "l1b", 'action', "present", 'file', strjoin(fullfile(l1b, l1bMembers), newline), 'source', url);
if mode == "never" || (isempty(missing) && mode ~= "force")
    if ~isempty(missing), rec.action = "missing (FetchData never)"; rec.file = ""; end
else
    if ~isfolder(l1b), mkdir(l1b); end
    tmp = [tempname '.tgz']; unpack = tempname;
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching GRACE-FO Level-1B %s (about 147 MB) from %s\n', l1bDay, url);
    try
        localDownload(tmp, url);
        mkdir(unpack);
        tarFile = gunzip(tmp, unpack);
        localTarMembers(tarFile{1}, l1bMembers, unpack);
        delete(tarFile{1});
        got = l1bMembers(isfile(fullfile(unpack, l1bMembers)));
        got = got(arrayfun(@(m) localHeadMatches(fullfile(unpack, m), '^header:'), got));
        for m = got, movefile(fullfile(unpack, m), fullfile(l1b, m), 'f'); end
        if numel(got) == numel(l1bMembers)
            rec.action = "fetched";
        else
            rec.action = "incomplete: the archive did not hold readable " + strjoin(setdiff(l1bMembers, got), ", ");
            rec.file = strjoin(fullfile(l1b, got), newline);
        end
    catch err
        rec.action = "failed: " + string(err.message);
        rec.file = "";
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, url);
    end
    if isfolder(unpack), rmdir(unpack, 's'); end
end
items(end+1) = rec;
d = struct('items', items, 'mode', mode);
end

function items = localFetchGroops(pool, base, files, mode)
%LOCALFETCHGROOPS The GROOPS Earth-radiation grids, into a dedicated DATAPOOL subfolder.
%   THE POOL IS NAMED, NEVER GUESSED. Without DATAPOOL this reports and does nothing: choosing a
%   directory for 2.5 GB on the user's behalf is not a service (R5), and the owner keeps the pool on
%   both the bridge and the CI runner.
%
%   NO SHA-256 PIN. Every other fetch here is pinned, and these are not: the files are 1.25 GB, the
%   publisher states plainly that it does not commit to keeping them updated, and a digest measured
%   once here would turn their next revision into a failure rather than a fact. What IS checked is
%   the GROOPS magic and type string, by reading the first bytes -- content, not a status code (F15).
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
if strlength(pool) == 0
    items(end+1) = struct('name', "groops", 'action', ...
        "skipped: DATAPOOL is not set -> set the DATAPOOL environment variable or pass 'DataPool'", ...
        'file', "", 'source', base);
    return
end
for k = 1:numel(files)
    folder = fullfile(pool, "EarthRadiationPressure");
    dst = fullfile(folder, files(k).name);
    url = base + "/" + files(k).folder + "/" + files(k).name;
    rec = struct('name', "groops " + extractBefore(files(k).name + ".", "."), 'action', "present", ...
        'file', string(dst), 'source', url);
    present = isfile(dst);
    if present && mode ~= "force"
        if ~localGroopsMagic(dst)
            rec.action = "failed: the file is not a GROOPS binary -> delete it and refetch";
            rec.file = "";
        end
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if mode == "never"
        rec.action = "missing (FetchGroops never)"; rec.file = "";
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if ~isfolder(folder), mkdir(folder); end
    tmp = [tempname '.dat'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching %s from %s\n', files(k).name, url);
    try
        localDownload(tmp, url);
        if localGroopsMagic(tmp)
            movefile(tmp, dst, 'f');
            rec.action = "fetched";
        else
            rec.action = "failed: the reply is not a GROOPS binary -> discarded";
            rec.file = "";
        end
    catch err
        rec.action = "failed: " + string(err.message);
        rec.file = "";
    end
    if rec.action ~= "fetched"
        fprintf('setupNonGravForces: %s -> fetch it by hand from %s\n', rec.action, url);
    end
    items(end+1) = rec; %#ok<AGROW>
    clear cleanup
end
end

function items = localFetchPool(pool, days, dstBase, f30Base, f30Files, itsgBase, itsgProducts, swarmHapi, swarmSats, mode)
%LOCALFETCHPOOL Observations of the reference days into DATAPOOL: Dst, F30, TU Graz's products, Swarm densities.
%   The pool is named, never guessed (localFetchGroops says why). Every reply is checked by its first
%   lines before it replaces anything -- for a .gz, those of what it decompresses to -- and a failure
%   is REPORTED per file, never raised. Kyoto keeps FINAL Dst only for months a few years old; a
%   recent day reports "failed" for its month rather than silently taking provisional values.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
if strlength(pool) == 0
    items(end+1) = struct('name', "pool", 'action', ...
        "skipped: DATAPOOL is not set -> set the DATAPOOL environment variable or pass 'DataPool'", ...
        'file', "", 'source', "");
    return
end
d = datetime(days, 'InputFormat', 'yyyy-MM-dd');
targets = struct('name', {}, 'url', {}, 'dst', {}, 'check', {});
% Dst: each day's month and the month before it, so a three-day history never stops at a month boundary
first = dateshift(d, 'start', 'month');
months = unique([first, first - calmonths(1)]);
for m = months
    ym = string(m, 'yyyyMM');
    yymm = extractAfter(ym, 2);
    targets(end+1) = struct('name', "dst " + ym, 'url', dstBase + "/" + ym + "/dst" + yymm + ".for.request", ...
        'dst', string(fullfile(pool, "SpaceWeather", "Dst", "dst" + yymm + ".for.request")), 'check', "^DST" + yymm); %#ok<AGROW>
end
for f = f30Files
    targets(end+1) = struct('name', "f30 " + extractBefore(f, "_observation"), 'url', f30Base + "/" + f, ...
        'dst', string(fullfile(pool, "SpaceWeather", "F30", f)), 'check', "^# File created by CLS"); %#ok<AGROW>
end
% TU Graz names GRACE-A/B "GRACE-1/-2" and GRACE-C/D "GRACEFO-1/-2"
for k = 1:numel(d)
    yyyy = string(year(d(k)));
    ymd = string(d(k), 'yyyy-MM-dd');
    if year(d(k)) < 2018, sats = ["GRACE-1" "GRACE-2"]; else, sats = ["GRACEFO-1" "GRACEFO-2"]; end
    for s = sats
        for p = itsgProducts
            file = s + "_" + p + "_" + ymd + ".txt.gz";
            targets(end+1) = struct('name', "itsg " + s + " " + p + " " + ymd, ...
                'url', itsgBase + "/" + s + "/" + p + "/" + yyyy + "/" + file, ...
                'dst', string(fullfile(pool, "ITSG", "operational", s, p, yyyy, file)), 'check', "^groops instrument"); %#ok<AGROW>
        end
    end
end
% Swarm DNSxPOD: the month of each day, from 2014 (Swarm was launched in November 2013)
targets = [targets, localSwarmTargets(pool, unique(first(year(first) >= 2014)), swarmHapi, swarmSats)];
items = localFetchTargets(items, targets, mode);
end

function t = localSwarmTargets(pool, months, swarmHapi, swarmSats)
%LOCALSWARMTARGETS One HAPI CSV per Swarm satellite and month (datetime first days of the months).
t = struct('name', {}, 'url', {}, 'dst', {}, 'check', {});
for m = months
    ym = string(m, 'yyyy-MM');
    span = "&start=" + string(m, 'yyyy-MM-dd') + "T00:00:00Z&stop=" + string(m + calmonths(1), 'yyyy-MM-dd') + "T00:00:00Z";
    for s = swarmSats
        ds = "SW_OPER_DNS" + s + "POD_2_";
        t(end+1) = struct('name', "swarm " + s + " " + ym, ...
            'url', swarmHapi + "?dataset=" + ds + span + "&format=csv&include=header", ...
            'dst', string(fullfile(pool, "Swarm", "DNSxPOD", string(year(m)), ds + ym + ".csv")), 'check', "^#\s*\{"); %#ok<AGROW>
    end
end
end

function items = localFetchTargets(items, targets, mode)
%LOCALFETCHTARGETS Fetch each target unless present, keep a reply only when its head matches.
for t = targets
    rec = struct('name', t.name, 'action', "present", 'file', t.dst, 'source', t.url);
    if isfile(t.dst) && mode ~= "force"
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if mode == "never"
        rec.action = "missing (FetchPool never)";
        rec.file = "";
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    folder = fileparts(t.dst);
    if ~isfolder(folder), mkdir(folder); end
    gz = endsWith(t.dst, ".gz");
    if gz, tmp = [tempname '.txt.gz']; else, tmp = [tempname '.txt']; end   % gunzip wants the extension
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching %s from %s\n', t.name, t.url);
    try
        if startsWith(t.url, "ftp://")
            localFtpGet(tmp, t.url);
        else
            localDownload(tmp, t.url);
        end
        if localPoolContent(tmp, t.check, gz)
            movefile(tmp, t.dst, 'f');
            rec.action = "fetched";
        else
            rec.action = "failed: the reply is not the expected product (its first lines do not match) -> discarded";
            rec.file = "";
        end
    catch err
        rec.action = "failed: " + string(err.message);
        rec.file = "";
    end
    if rec.action ~= "fetched"
        fprintf('setupNonGravForces: %s -> fetch it by hand from %s\n', rec.action, t.url);
    end
    items(end+1) = rec; %#ok<AGROW>
    clear cleanup
end
end

function items = localFetchPoolMonths(pool, months, l1bBase, products, swarmHapi, swarmSats, eopUrl, swUrl, mode)
%LOCALFETCHPOOLMONTHS Whole months into DATAPOOL: GRACE-C/D Level-1B members, Swarm, SW-All and finals2000A.
%   Everything a month run of tools/dragScaleMonth.m reads. The two complete files are judged by
%   COVERAGE, not existence: one that ends before the last month's last day is refetched, since an
%   old copy is a file that exists and is wrong. A Level-1B day is present when all its members are;
%   otherwise its archive (about 150 MB) is fetched, the members are cut out with localTarMembers and
%   kept only if each starts with an L1B 'header:', and the archive is deleted.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
if isempty(months), return, end
if strlength(pool) == 0
    items(end+1) = struct('name', "pool months", 'action', ...
        "skipped: DATAPOOL is not set -> set the DATAPOOL environment variable or pass 'DataPool'", 'file', "", 'source', "");
    return
end
m = unique(datetime(months, 'InputFormat', 'yyyy-MM'));
lastDay = dateshift(max(m), 'end', 'month');
% complete space weather and EOP, refetched when stale
whole = struct('name', {"pool spaceweather", "pool eop"}, 'url', {swUrl, eopUrl}, ...
    'dst', {string(fullfile(pool, "SpaceWeather", "SW-All.csv")), string(fullfile(pool, "EOP", "finals2000A.all"))}, ...
    'check', {'^DATE,BSRN,ND,KP1', '^\s*\d{1,2}\s*\d{1,2}\s*\d{1,2}\s+\d{5}\.\d{2}\s+[IP]'});
kinds = ["sw" "eop"];
for k = 1:numel(whole)
    use = mode;
    if mode == "ifMissing" && isfile(whole(k).dst) && ~(localLastDate(whole(k).dst, kinds(k)) >= lastDay)
        use = "force";
    end
    items = localFetchTargets(items, whole(k), use);
end
items = localFetchTargets(items, localSwarmTargets(pool, m(year(m) >= 2014), swarmHapi, swarmSats), mode);
% GRACE-FO Level-1B, both satellites, every day of every month
for mm = m
    for day = mm:caldays(1):dateshift(mm, 'end', 'month')
        ymd = string(day, 'yyyy-MM-dd');
        dst = fullfile(pool, "GRACEFO", "L1B", "RL04", string(year(day)), ymd);
        members = reshape(products' + "_" + ymd + "_" + ["C" "D"] + "_04.txt", 1, []);
        url = l1bBase + "/" + string(year(day)) + "/gracefo_1B_" + ymd + "_RL04.ascii.noLRI.tgz";
        rec = struct('name', "pool l1b " + ymd, 'action', "present", 'file', string(dst), 'source', url);
        missing = members(~isfile(fullfile(dst, members)));
        if isempty(missing) && mode ~= "force"
            items(end+1) = rec; %#ok<AGROW>
            continue
        end
        if mode == "never"
            rec.action = "missing (FetchPool never)"; rec.file = "";
            items(end+1) = rec; %#ok<AGROW>
            continue
        end
        tmp = [tempname '.tgz']; unpack = tempname;
        cleanup = onCleanup(@() localQuietDelete(tmp));
        fprintf('setupNonGravForces: fetching GRACE-FO Level-1B %s (about 150 MB) from %s\n', ymd, url);
        try
            localDownload(tmp, url);
            mkdir(unpack);
            tarFile = gunzip(tmp, unpack);
            localTarMembers(tarFile{1}, members, unpack);
            delete(tarFile{1});
            got = members(isfile(fullfile(unpack, members)));
            got = got(arrayfun(@(f) localHeadMatches(fullfile(unpack, f), '^header:'), got));
            if ~isfolder(dst), mkdir(dst); end
            for f = got, movefile(fullfile(unpack, f), fullfile(dst, f), 'f'); end
            if numel(got) == numel(members)
                rec.action = "fetched";
            else
                rec.action = "incomplete: the archive did not hold readable " + strjoin(setdiff(members, got), ", ");
            end
        catch err
            rec.action = "failed: " + string(err.message);
            rec.file = "";
            fprintf('setupNonGravForces: %s -> fetch it by hand from %s\n', rec.action, url);
        end
        if isfolder(unpack), rmdir(unpack, 's'); end
        localQuietDelete(tmp);
        clear cleanup
        items(end+1) = rec; %#ok<AGROW>
    end
end
end

function t = localLastDate(file, kind)
%LOCALLASTDATE The date of a complete file's last record: SW-All "yyyy-mm-dd,...", finals2000A "yymmdd".
t = NaT;
fid = fopen(file, 'r');
if fid < 0, return, end
closer = onCleanup(@() fclose(fid));
fseek(fid, 0, 'eof');
n = ftell(fid);
fseek(fid, max(0, n - 8192), 'bof');
lines = splitlines(string(fread(fid, [1 Inf], '*char')));
lines = lines(strlength(strtrim(lines)) > 0);
if isempty(lines), return, end
last = lines(end);
if kind == "sw"
    tok = regexp(last, '^(\d{4}-\d{2}-\d{2}),', 'tokens', 'once');
    if ~isempty(tok), t = datetime(tok{1}, 'InputFormat', 'yyyy-MM-dd'); end
else
    tok = regexp(last, '^\s*(\d{1,2})\s*(\d{1,2})\s*(\d{1,2})\s', 'tokens', 'once');
    if ~isempty(tok)
        y = str2double(tok{1}); y = y + 1900 + 100 * (y < 73);
        t = datetime(y, str2double(tok{2}), str2double(tok{3}));
    end
end
end

function localMustBeMonths(m)
%LOCALMUSTBEMONTHS PoolMonths are "yyyy-mm" from 2018-06, when GRACE-FO's Level-1B begins.
if isempty(m), return, end
if ~all(matches(m, digitsPattern(4) + "-" + digitsPattern(2)))
    error('setupNonGravForces:poolMonths', 'PoolMonths must be "yyyy-mm" strings -> e.g. "2024-01"');
end
if any(datetime(m, 'InputFormat', 'yyyy-MM') < datetime(2018, 6, 1))
    error('setupNonGravForces:poolMonths', 'PoolMonths start at 2018-06 (GRACE-FO Level-1B) -> pass a later month');
end
end

function localFtpGet(tmp, url)
%LOCALFTPGET One file by anonymous FTP: websave speaks HTTP only, and CLS publishes F30 by FTP only.
tok = regexp(char(url), '^ftp://([^/]+)(/.*)/([^/]+)$', 'tokens', 'once');
if isempty(tok)
    error('setupNonGravForces:ftpUrl', 'not an ftp://host/path/file URL: %s', url);
end
here = tempname;
mkdir(here);
gone = onCleanup(@() rmdir(here, 's'));
f = ftp(tok{1});
closer = onCleanup(@() close(f));
cd(f, tok{2});
mget(f, tok{3}, here);
clear closer
movefile(fullfile(here, tok{3}), tmp, 'f');
clear gone
end

function tf = localPoolContent(file, pattern, gz)
%LOCALPOOLCONTENT The first lines of a fetched observation -- decompressed first when it is a .gz -- match.
if ~gz
    tf = localHeadMatches(file, pattern);
    return
end
tf = false;
here = tempname;
mkdir(here);
gone = onCleanup(@() rmdir(here, 's'));
try
    out = gunzip(file, here);
catch
    return          % not gzip: an HTML error page, typically
end
tf = ~isempty(out) && localHeadMatches(out{1}, pattern);
clear gone
end

function tf = localGroopsMagic(file)
%LOCALGROOPSMAGIC The first bytes are 'B' plus a length-prefixed type name starting "gridded".
%   A status code is not a file (CLAUDE.md), and an HTML error page saved as .dat is exactly the
%   failure this catches.
tf = false;
fid = fopen(file, 'r');
if fid < 0, return; end
closer = onCleanup(@() fclose(fid));
if fread(fid, 1, '*char') ~= 'B', return; end
n = fread(fid, 1, 'uint64');
if isempty(n) || n < 7 || n > 64, return; end
tf = startsWith(string(fread(fid, n, '*char')'), "gridded");
end

function items = localFetchKernel(root, file, sha256, mode)
%LOCALFETCHKERNEL The DE440 planetary ephemeris, fetched through referenceSystem's catalogue.
%   TWO JOBS, ONE AUTHORITY EACH (R1). referenceSystem.downloadKernels owns where the kernel comes
%   from -- it says so in its own help and this toolbox does not repeat the URL. What belongs here
%   is the pin: the SHA-256 of the file this project has actually tested against, so a silently
%   reissued kernel is a finding and not a mystery.
%
%   It is fetched and never committed: 32.7 MB of binary that NAIF publishes and keeps current is
%   the same case as the MSIS and HWM14 sources (O-34), and the same answer.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
folder = fullfile(root, 'data', 'kernels');
dst = fullfile(folder, file);
rec = struct('name', "de440", 'action', "present", 'file', string(dst), ...
    'source', "referenceSystem.downloadKernels(Set=minimal, Span=short)");
present = isfile(dst);
if present && mode ~= "force"
    % CONTENT, NOT EXISTENCE. A pin checked only at download time protects the first run on a
    % fresh machine and nothing after it, which is every run on every developer's host. Verified
    % 2026-09-16 by fault injection: with the digest checked only on download, a deliberately wrong
    % pin changed nothing and the gate stayed green. Hashing 33 MB costs about 0.2 s.
    try
        got = localSha256(dst);
    catch err
        got = "";
        rec.action = "failed: " + string(err.message);
        rec.file = "";
    end
    if strlength(got) > 0 && got ~= sha256
        rec.action = "failed: the file on disk has SHA-256 " + got + ", not the pinned " + sha256 + ...
            " -> it was NOT used; delete it and refetch, or verify it by hand before repinning";
        rec.file = "";
    end
    items(end+1) = rec;
    return
end
if mode == "never"
    rec.action = "missing (FetchKernels never)"; rec.file = "";
    items(end+1) = rec;
    return
end
if strlength(string(which('referenceSystem.downloadKernels'))) == 0
    % The installer does not require referenceSystem to be on the path yet; say so instead of
    % failing, because every other fetch here works without it.
    rec.action = "skipped: referenceSystem is not on the path -> run setupReferenceSystem first";
    rec.file = "";
    items(end+1) = rec;
    return
end
if ~isfolder(folder), mkdir(folder); end
fprintf('setupNonGravForces: fetching %s through referenceSystem.downloadKernels (about 33 MB)\n', file);
try
    % READ THE MANIFEST BEFORE THE FILE. downloadKernels does not throw on a failed fetch: it warns
    % and returns a "FAILED | reason" row, and writes nothing. On CI run 38733381139 NAIF closed the
    % connection, the file was never produced, and hashing it reported "Invalid file identifier" in
    % place of the network error -- the gate went red for the right reason under the wrong name.
    manifest = referenceSystem.downloadKernels('Set', "minimal", 'Span', "short", 'Mode', "download", ...
        'CacheDir', string(folder), 'Overwrite', mode == "force");
    row = manifest(manifest.File == string(file), :);
    if ~isempty(row) && startsWith(row.Status(1), "FAILED")
        error('setupNonGravForces:download', '%s', ...
            "download: " + strtrim(extractAfter(row.Status(1), "FAILED |")) + ...
            " -> nothing was written; rerun setupNonGravForces('FetchKernels', ""force"") when NAIF answers");
    end
    if ~isfile(dst)
        error('setupNonGravForces:download', '%s', "download: referenceSystem.downloadKernels reported " + ...
            "status """ + strjoin(string(row.Status), ", ") + """ but wrote no " + string(file) + " in " + string(folder) + ...
            " -> fetch it by hand into that folder, or rerun when NAIF answers");
    end
    got = localSha256(dst);
    if got == sha256
        rec.action = "fetched";
    else
        % DELETE IT. downloadKernels writes straight to the destination, so a file that fails its
        % check is sitting exactly where every consumer looks; on 2026-09-16 the tests read one and
        % passed, and the run reported "failed" and GREEN in the same log. An unverified file must
        % not survive the verification that rejected it.
        localQuietDelete(dst);
        rec.action = "failed: SHA-256 " + got + " is not the pinned " + sha256 + ...
            " -> the file was DELETED; NAIF may have reissued the kernel, so verify it by hand before repinning";
        rec.file = "";
    end
catch err
    rec.action = "failed: " + string(err.message);
    rec.file = "";
end
if rec.action ~= "fetched" && rec.action ~= "present"
    fprintf('setupNonGravForces: %s\n', rec.action);
end
items(end+1) = rec;
end

function items = localFetchCeres(root, base, files, sha256, mode)
%LOCALFETCHCERES The CERES EBAF files, each kept only if its SHA-256 is ESGF's.
items = struct('name', {}, 'action', {}, 'file', {}, 'source', {});
folder = fullfile(root, 'data', 'ceres');
for k = 1:numel(files)
    variable = extractBefore(files(k), "_");
    url = base + "/" + variable + "/gn/v20260220/" + files(k);
    dst = fullfile(folder, files(k));
    rec = struct('name', "ceres " + variable, 'action', "present", 'file', string(dst), 'source', url);
    present = isfile(dst);
    if mode == "never" || (present && mode ~= "force")
        if ~present, rec.action = "missing (FetchData never)"; rec.file = ""; end
        items(end+1) = rec; %#ok<AGROW>
        continue
    end
    if ~isfolder(folder), mkdir(folder); end
    tmp = [tempname '.nc'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    fprintf('setupNonGravForces: fetching CERES EBAF %s from %s\n', variable, url);
    try
        localDownload(tmp, url);
        got = localSha256(tmp);
        if got == sha256(k)
            movefile(tmp, dst, 'f');
            rec.action = "fetched";
        else
            rec.action = "failed: SHA-256 " + got + " is not ESGF's " + sha256(k) + " -> the reply was discarded";
            rec.file = "";
        end
    catch err
        rec.action = "failed: " + string(err.message);
        rec.file = "";
    end
    if rec.action ~= "fetched"
        fprintf('setupNonGravForces: %s -> download it by hand from %s\n', rec.action, url);
    end
    items(end+1) = rec; %#ok<AGROW>
    clear cleanup
end
end

function localTarMembers(tarFile, names, outDir)
%LOCALTARMEMBERS Write the members of a POSIX (ustar or GNU) tar whose base names are in names.
% A tar is a sequence of 512-byte headers, each followed by its data padded to 512 bytes; an
% all-zero header ends it. GNU long names arrive as a preceding 'L' entry. Only regular files
% ('0' or NUL) are written, under their base name, so a folder inside the archive does not matter.
fid = fopen(tarFile, 'r');
if fid < 0
    error('setupNonGravForces:tar', 'cannot open %s', tarFile);
end
closer = onCleanup(@() fclose(fid));
longName = "";
while true
    h = fread(fid, 512, '*uint8')';
    if numel(h) < 512 || all(h == 0), break; end
    name = string(deblank(char(h(1:100))));
    prefix = string(deblank(char(h(346:500))));
    if strlength(prefix) > 0, name = prefix + "/" + name; end
    digits = char(h(125:136));
    digits = digits(digits >= '0' & digits <= '7');       % octal, NUL- or blank-terminated, sometimes blank-padded
    sz = 0;
    if ~isempty(digits), sz = base2dec(digits, 8); end
    type = char(h(157));
    blocks = ceil(sz / 512) * 512;
    if type == 'L'
        data = fread(fid, blocks, '*uint8')';
        longName = string(deblank(char(data(1:sz))));
        continue
    end
    if strlength(longName) > 0, name = longName; longName = ""; end
    [~, base, ext] = fileparts(name);
    if (type == '0' || type == char(0)) && any(base + ext == names)
        data = fread(fid, sz, '*uint8');
        out = fopen(fullfile(outDir, base + ext), 'w');
        fwrite(out, data, 'uint8');
        fclose(out);
        fseek(fid, blocks - sz, 'cof');
    else
        fseek(fid, blocks, 'cof');
    end
end
end

function tf = localHeadMatches(file, pattern)
%LOCALHEADMATCHES The first ten lines of a file contain the pattern: the content check a status code is not.
tf = false;
fid = fopen(file, 'r');
if fid < 0, return; end
closer = onCleanup(@() fclose(fid));
for k = 1:10
    line = fgetl(fid);
    if ~ischar(line), return; end
    if ~isempty(regexp(line, pattern, 'once'))
        tf = true;
        return
    end
end
end

function localQuietDelete(f)
%LOCALQUIETDELETE Remove a temporary file without turning a cleanup into an error.
if isfile(f)
    try, delete(f); catch, end
end
end

function rec = localFetchArchive(name, dst, marker, url, sha256, mode, notice)
%LOCALFETCHARCHIVE A pinned tarball of model source, unpacked into dst after its digest matched.
rec = struct('name', name, 'action', "present", 'dir', string(dst), 'source', url, 'message', "");
here = isfile(fullfile(dst, marker));
if mode == "never" || (here && mode ~= "force")
    if ~here, rec.action = "skipped"; end
    return
end
fprintf('setupNonGravForces: fetching %s from %s\n', name, url);
tmp = [tempname '.tar.gz'];
cleanup = onCleanup(@() localQuietDelete(tmp));
try
    localDownload(tmp, url);
catch err
    rec.action = "failed"; rec.message = "download failed: " + string(err.message) + " -> fetch " + url + " by hand and unpack it into " + dst;
    fprintf('setupNonGravForces: %s\n', rec.message);
    return
end
got = localSha256(tmp);
if got ~= sha256
    rec.action = "failed";
    rec.message = "the archive's SHA-256 is " + got + ", the pin is " + sha256 + " -> the release changed or the download is not the archive; nothing was unpacked";
    fprintf('setupNonGravForces: %s\n', rec.message);
    return
end
if isfolder(dst), rmdir(dst, 's'); end
untar(tmp, dst);
fprintf('setupNonGravForces: %s\n', notice);
rec.action = "fetched";
rec.message = notice;
end

function rec = localFetchFiles(name, dst, files, mode)
%LOCALFETCHFILES Model source files fetched one by one into dst, each kept only if its SHA-256 is the pin.
rec = struct('name', name, 'action', "present", 'dir', string(dst), 'source', strjoin([files.url], " "), 'message', "");
here = all(arrayfun(@(f) isfile(fullfile(dst, f.name)), files));
if mode == "never" || (here && mode ~= "force")
    if ~here, rec.action = "skipped"; end
    return
end
if ~isfolder(dst), mkdir(dst); end
fprintf('setupNonGravForces: fetching %s (%d files)\n', name, numel(files));
for f = files
    tmp = [tempname '.part'];
    cleanup = onCleanup(@() localQuietDelete(tmp));
    try
        localDownload(tmp, f.url);
    catch err
        rec.action = "failed"; rec.message = f.name + ": " + string(err.message) + " -> fetch " + f.url + " by hand into " + dst;
        fprintf('setupNonGravForces: %s\n', rec.message);
        return
    end
    got = localSha256(tmp);
    if got ~= f.sha256
        rec.action = "failed"; rec.message = f.name + ": SHA-256 " + got + " is not the pin " + f.sha256 + " -> the reply was discarded";
        fprintf('setupNonGravForces: %s\n', rec.message);
        return
    end
    movefile(tmp, fullfile(dst, f.name), 'f');
    clear cleanup
end
rec.action = "fetched";
end

function localDownload(tmp, url)
%LOCALDOWNLOAD Every fetch of this file: websave with a 60 s timeout, retried once after 5 s.
%   MATLAB's default weboptions timeout is 5 s, and with it CI run 38786016017 lost NRLMSIS 2.1 to
%   "timed out after 5.000 seconds" from map.nrl.navy.mil; earlier that day NAIF closed a transfer
%   mid-file (run 38733381139). 60 s is a slow research server's first byte, not a transfer limit:
%   the timeout bounds the wait for data, and a 147 MB archive keeps data flowing. One retry covers
%   a dropped connection or a 503; a server that fails twice is down, and the caller REPORTS it with
%   both causes. Nothing here checks content: every caller verifies what arrives (SHA-256 pins, the
%   product's first lines, the GROOPS magic) after this returns, retry or not.
opts = weboptions('Timeout', 60);
first = "";
for attempt = 1:2
    try
        websave(tmp, char(url), opts);
        return
    catch err
        if attempt == 1
            first = string(err.message);
            pause(5);
        else
            error('setupNonGravForces:download', '%s', "two attempts failed: (1) " + first + " (2) " + string(err.message));
        end
    end
end
end

function h = localSha256(file)
md = java.security.MessageDigest.getInstance('SHA-256');
[fid, why] = fopen(file, 'r');
if fid < 0
    % Refuse here, BEFORE registering fclose: an onCleanup on fid -1 raises "Invalid file
    % identifier", which is what reached the gate in place of the real cause (CI run 38733381139).
    error('setupNonGravForces:sha256', 'cannot open %s to hash it (%s) -> the file is missing or unreadable', ...
        file, why);
end
closer = onCleanup(@() fclose(fid));
md.update(fread(fid, Inf, '*uint8'));
h = string(lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, [])));
end

function m = localFetchMsis(root, url, commit, mode)
%LOCALFETCHMSIS Obtain the NRLMSISE-00 C source if it is not already here.
% Absence is reported rather than raised (R4/R5): buildMex is the unit that refuses, and one
% cause should produce one error and not two.
m = struct('dir', "", 'action', "skipped", 'url', url, 'commit', commit, 'message', "");
env = string(getenv('NGF_NRLMSISE00'));
if strlength(env) > 0
    m.dir = env;
    m.action = "skipped";
    m.message = "NGF_NRLMSISE00 points at a copy, so nothing is fetched";
    return
end
dst = string(fullfile(root, 'mex', 'nrlmsise00'));
m.dir = dst;
here = isfile(fullfile(dst, 'nrlmsise-00.c'));
if mode == "never" || (here && mode ~= "force")
    if here, m.action = "present"; else, m.action = "skipped"; end
    return
end
if here && mode == "force"
    rmdir(dst, 's');
end
% One line, whatever Verbose says: a network action has to be visible.
fprintf('setupNonGravForces: fetching NRLMSISE-00 from %s (commit %s)\n', url, extractBefore(commit, 8));
% No --depth: this server speaks dumb HTTP, which does not support shallow clones (measured
% 2026-09-09). The repository is a few hundred kilobytes, so a full clone costs nothing.
[st, out] = system(sprintf('git clone --quiet "%s" "%s"', url, dst));
if st == 0
    [st, out] = system(sprintf('git -C "%s" checkout --quiet %s', dst, commit));
end
if st ~= 0
    if isfolder(dst), rmdir(dst, 's'); end        % never leave a half-clone for buildMex to digest
    m.action = "failed";
    m.message = "git failed: " + string(strtrim(out)) + " -> clone it by hand: git clone " + url + " mex/nrlmsise00";
    fprintf('setupNonGravForces: %s\n', m.message);
    return
end
m.action = "fetched";
end
