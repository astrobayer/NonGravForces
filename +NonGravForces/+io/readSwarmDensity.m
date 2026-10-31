function [epoch, lla, rho, aux] = readSwarmDensity(files, options)
%READSWARMDENSITY Thermosphere density from VirES HAPI CSV: Swarm POD-only, CHAMP/GRACE/GRACE-FO accelerometer.
%
%   [epoch, lla, rho, aux] = NonGravForces.io.readSwarmDensity(file)
%   [epoch, lla, rho, aux] = NonGravForces.io.readSwarmDensity(files, 'TimeRange', [mjd0 mjd1])
%
%   Reads ESA's Level-2 thermosphere density of Swarm A, B or C derived from
%   precise orbit determination ALONE (no accelerometer), as the VirES for
%   Swarm HAPI server serves it: CSV with the HAPI JSON header in '#' lines,
%   which is what setupNonGravForces('FetchPool') stores in DATAPOOL and what
%   tools/cutFixtures.py cuts into data/samples/swarm. The header is decoded
%   and the columns are found BY NAME, and their units are checked against
%   what this reader converts from; a reordered or relabelled product is
%   refused rather than read in the wrong column.
%
%   The density is an observation independent of GRACE-FO's accelerometer,
%   which is what makes it usable as a density SCALE for GRACE-FO's drag
%   (HANDOVER O-19): rho_Swarm / rho_model at Swarm's positions.
%
%   The same server serves TU Delft's ACCELEROMETER-derived densities of CHAMP,
%   GRACE and GRACE-FO (ESA TOLEOS, datasets CH_OPER_DNS_ACC_2_, GR_OPER_DNS<1|2>ACC_2_,
%   GF_OPER_DNS<1|2>ACC_2_; CC BY 4.0) in exactly this layout -- the same nine
%   parameters and units, at 10 s -- and this reader reads them too (2026-10-28,
%   O-19 step 3: densities of solar cycle 23 for tools/densityDriverFit.m). The name
%   stays: renaming a public unit breaks its callers. TOLEOS marks a missing value
%   with the fill 9.99e+32 although its header declares no fill; any value >= 1e30
%   is returned as NaN, and a record whose density is NaN is not valid.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more HAPI CSV files of ONE dataset (SW_OPER_DNS<A|B|C>POD_2_,
%                                 CH_OPER_DNS_ACC_2_, GR_OPER_DNS<1|2>ACC_2_, GF_OPER_DNS<1|2>ACC_2_);
%                                 concatenated and sorted by time, duplicates of an epoch kept once.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type    Units     Description
%   ----        ----  ----    -----     -----------
%   TimeRange   1x2   double  MJD(UTC)  Keep records with mjd0 <= t <= mjd1. Default [] (all).
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units      Description
%   ----     ----  ----    -----      -----------
%   epoch    1x1   struct  -          referenceSystem epoch struct, K entries, scale "UTC"
%   lla      3xK   double  rad,rad,m  Geodetic latitude, longitude, height (GRS80, as the
%                                     product states; see ACCURACY for WGS-84)
%   rho      1xK   double  kg/m^3     Thermosphere mass density (NaN where the product has its fill)
%   aux      1x1   struct  -          .rhoOrbitMean (1xK kg/m^3, the product's orbit
%                                     average, NaN at the fill), .lst_h (1xK local solar time, h),
%                                     .valid (1xK logical, validity_flag == 0 and the density not
%                                     the fill), .satellite ("A"|"B"|"C" for Swarm, "CHAMP",
%                                     "GRACE-A"|"GRACE-B", "GRACE-C"|"GRACE-D" for GRACE-FO 1|2),
%                                     .dataset, .header (decoded HAPI header of the first
%                                     file), .files, .numRecords, .provenance
%
% ACCURACY
%   Exact transcription: degrees to radians by deg2rad, the time tag's milliseconds kept,
%   the calendar date to MJD by Fliegel & Van Flandern's integer formula (exact 1901-2099).
%   The positions are GRS80 geodetic and env.density takes WGS-84 geodetic: the two
%   ellipsoids differ by 0.1 mm in semi-minor axis, so the heights agree to 1e-4 m and
%   the latitudes to 1e-9 rad, far below anything a density depends on. The density
%   itself is ESA's: derived from POD-only non-gravitational accelerations with a
%   satellite aerodynamic model, nominal 30 s sampling, effective resolution of the
%   order of 20 min, uncertainty stated as 30 % of the model variance (Swarm L2 product
%   specification); it is NOT smoothed or filtered here. TOLEOS's densities come from the
%   calibrated accelerometers with TU Delft's high-fidelity geometry models, version 02
%   (March et al. 2019, 2021), and are likewise passed through unchanged.
%
% ERRORS
%   NonGravForces:readSwarmDensity:inputRange    TimeRange is not [mjd0 mjd1] with mjd0 <= mjd1
%   NonGravForces:readSwarmDensity:dataFormat     no HAPI JSON header, an undecodable one, or a
%                                                 status other than 1200 (OK)
%   NonGravForces:readSwarmDensity:dataColumn     the first parameter is not Timestamp, a required
%                                                 parameter is missing from the header, or a data
%                                                 line does not have the header's columns
%   NonGravForces:readSwarmDensity:dataUnits      a parameter's units are not the ones converted here
%   NonGravForces:readSwarmDensity:dataSatellite  the files are of different datasets, or of a
%                                                 dataset this reader does not know
%   NonGravForces:readSwarmDensity:dataEmpty      no record survives TimeRange
%
% REFERENCES
%   ESA Swarm Level 2 Product Specification, SW_DNSxPOD_2_ (density, density_orbitmean,
%   validity_flag); van den IJssel, Encarnacao, Doornbos, Visser (2015), Adv. Space Res. 56,
%   1015-1029 (thermosphere densities from Swarm GPS); March, Visser, Doornbos et al. (2019),
%   Adv. Space Res. 63, 213-238 (CHAMP, GOCE, GRACE, Swarm densities with high-fidelity
%   geometries); HAPI 3.0 specification (header and CSV).
%
% EXAMPLE
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'swarm', ...
%       'SW_OPER_DNSBPOD_2_2020-01-01_2020-01-03.csv');
%   [ep, lla, rho, aux] = NonGravForces.io.readSwarmDensity(f); median(rho(aux.valid))
%
% See also NonGravForces.env.density, setupNonGravForces, NonGravForces.io.spaceWeather.
%
% KEYWORDS Swarm CHAMP GRACE TOLEOS density POD accelerometer thermosphere HAPI VirES drag scale
% Created 2026-10-27 by Claude (Anthropic) for Matthias Weigelt.
arguments
    files (1,:) string {mustBeNonempty, mustBeFile}
    options.TimeRange double = []
end
unit = "readSwarmDensity";
if ~isempty(options.TimeRange) && (numel(options.TimeRange) ~= 2 || options.TimeRange(2) < options.TimeRange(1))
    NonGravForces.fail(unit, "inputRange", "TimeRange must be [mjd0 mjd1] with mjd0 <= mjd1 -> pass two UTC MJDs");
end
% name -> the units this reader converts from
need = ["Latitude_GD" "Longitude_GD" "Height_GD" "local_solar_time" "density" "density_orbitmean" "validity_flag"];
units = ["deg" "deg" "m" "h" "kg/m3" "kg/m3" ""];
mjd = zeros(0, 1); sec = zeros(0, 1); val = zeros(0, numel(need));
dataset = ""; header = [];
for f = files
    [hdr, body] = localSplit(unit, f);
    if dataset == ""
        dataset = string(hdr.x_dataset); header = hdr;
        satellite = localSatellite(unit, f, dataset);
    elseif string(hdr.x_dataset) ~= dataset
        NonGravForces.fail(unit, "dataSatellite", "%s is %s, the first file %s -> read one satellite per call", f, hdr.x_dataset, dataset);
    end
    par = hdr.parameters;                         % a cell when the parameters carry different fields
    if ~iscell(par), par = num2cell(par); end
    names = string(cellfun(@(p) p.name, par, 'UniformOutput', false));
    if names(1) ~= "Timestamp"
        NonGravForces.fail(unit, "dataColumn", "%s: the first parameter is %s, not Timestamp -> this is not a VirES HAPI density file", f, names(1));
    end
    col = zeros(1, numel(need));
    for k = 1:numel(need)
        j = find(names == need(k), 1);
        if isempty(j)
            NonGravForces.fail(unit, "dataColumn", "%s: parameter %s is missing from the header -> the product changed; read it before adapting this reader", f, need(k));
        end
        u = par{j}.units;
        if isempty(u), u = ""; end
        if string(u) ~= units(k)
            NonGravForces.fail(unit, "dataUnits", "%s: %s is in '%s', this reader converts from '%s' -> the product changed; adapt the conversion", f, need(k), string(u), units(k));
        end
        col(k) = j;
    end
    nPar = numel(names);
    fmt = ['%4f-%2f-%2fT%2f:%2f:%fZ' repmat(',%f', 1, nPar - 1)];
    C = textscan(body, fmt, 'CollectOutput', true);
    M = C{1};
    if isempty(M) || size(M, 2) ~= 6 + nPar - 1 || any(isnan(M(:, 1)))
        NonGravForces.fail(unit, "dataColumn", "%s: the data lines do not have the header's %d columns -> check the file", f, nPar);
    end
    nLines = numel(strfind(body, newline)) + 1;  % localSplit joined the non-empty data lines with newline
    if nLines ~= size(M, 1)
        NonGravForces.fail(unit, "dataColumn", "%s: %d data lines but %d parsed -> a line has a different number of columns", f, nLines, size(M, 1));
    end
    mjd = [mjd; localMjd(M(:, 1), M(:, 2), M(:, 3))]; %#ok<AGROW>
    sec = [sec; M(:, 4) * 3600 + M(:, 5) * 60 + M(:, 6)]; %#ok<AGROW>
    val = [val; M(:, 6 + col - 1)]; %#ok<AGROW> parameter j is column 6 + j - 1 (the time takes six)
end
[~, order] = unique([mjd sec], 'rows');          % sorted, an epoch present in two files kept once
mjd = mjd(order); sec = sec(order); val = val(order, :);
if ~isempty(options.TimeRange)
    t = mjd + sec / 86400;
    keep = t >= options.TimeRange(1) & t <= options.TimeRange(2);
    mjd = mjd(keep); sec = sec(keep); val = val(keep, :);
end
if isempty(mjd)
    NonGravForces.fail(unit, "dataEmpty", "no record of %s in the files or in TimeRange [%s] -> widen TimeRange or read the month that holds it", ...
        dataset, num2str(options.TimeRange));
end
epoch = referenceSystem.daySecondToEpoch([mjd sec], "UTC");
lla = [deg2rad(val(:, 1))'; deg2rad(val(:, 2))'; val(:, 3)'];
r = val(:, 5:6); r(r >= 1e30) = NaN; val(:, 5:6) = r;   % TOLEOS's fill 9.99e+32 (its header declares none)
rho = val(:, 5)';
aux.rhoOrbitMean = val(:, 6)';
aux.lst_h = val(:, 4)';
aux.valid = val(:, 7)' == 0 & ~isnan(rho);
aux.satellite = satellite;
aux.dataset = dataset;
aux.header = header;
aux.files = files;
aux.numRecords = numel(mjd);
aux.provenance = NonGravForces.provenance("io.readSwarmDensity", 'Files', files, 'Options', options);
end

function s = localSatellite(unit, file, dataset)
%LOCALSATELLITE The spacecraft of a known dataset: Swarm A/B/C, CHAMP, GRACE-A/B, GRACE-C/D (GRACE-FO 1/2).
t = regexp(dataset, '^(SW_OPER_DNS([ABC])POD_2_|CH_OPER_DNS_ACC_2_|GR_OPER_DNS([12])ACC_2_|GF_OPER_DNS([12])ACC_2_)$', 'tokens', 'once');
if isempty(t)
    NonGravForces.fail(unit, "dataSatellite", "%s is dataset %s, which this reader does not know -> read Swarm DNSxPOD or TOLEOS DNSxACC files", file, dataset);
end
d = char(dataset);
switch d(1:2)
    case 'SW', s = string(d(12));
    case 'CH', s = "CHAMP";
    case 'GR', s = "GRACE-" + string(char('A' + (d(12) - '1')));
    otherwise, s = "GRACE-" + string(char('C' + (d(12) - '1')));
end
end

function mjd = localMjd(y, m, d)
%LOCALMJD Integer MJD of a Gregorian calendar date, 1901-2099 (Fliegel & Van Flandern 1968, shifted to MJD).
mjd = 367 * y - floor(7 * (y + floor((m + 9) / 12)) / 4) + floor(275 * m / 9) + d - 678987;
end

function [hdr, body] = localSplit(unit, file)
%LOCALSPLIT The decoded HAPI header ('#' lines) and the data lines of one file.
txt = fileread(file);
lines = splitlines(string(txt));
isHead = startsWith(lines, "#");
if ~any(isHead) || ~startsWith(strtrim(extractAfter(lines(find(isHead, 1)), 1)), "{")
    NonGravForces.fail(unit, "dataFormat", "%s has no HAPI JSON header ('#{' first) -> fetch it through setupNonGravForces('FetchPool')", file);
end
try
    hdr = jsondecode(char(strjoin(extractAfter(lines(isHead), 1), newline)));
catch err
    NonGravForces.fail(unit, "dataFormat", "%s: the HAPI header is not JSON (%s) -> the file is truncated or not a HAPI reply", file, err.message);
end
if ~isfield(hdr, 'status') || hdr.status.code ~= 1200 || ~isfield(hdr, 'parameters') || ~isfield(hdr, 'x_dataset')
    NonGravForces.fail(unit, "dataFormat", "%s: the HAPI header has no OK status (1200), parameters or x_dataset -> refetch it", file);
end
data = lines(~isHead & strlength(strtrim(lines)) > 0);
body = char(strjoin(data, newline));
end
