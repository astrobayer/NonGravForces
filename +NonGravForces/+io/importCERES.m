function [cer, aux] = importCERES(files, options)
%IMPORTCERES CERES EBAF TOA monthly fluxes (obs4MIPs netCDF) as the grids env.earthRadiance("ceres") evaluates.
%
%   cer = NonGravForces.io.importCERES(files)
%   cer = NonGravForces.io.importCERES(files, 'Period', [datetime(2020,1,1) datetime(2020,1,31)])
%   [cer, aux] = NonGravForces.io.importCERES(___)
%
%   Reads the three obs4MIPs files of the CERES Energy Balanced and Filled TOA
%   product -- outgoing shortwave (rsut), outgoing longwave (rlut) and incident
%   shortwave (rsdt), monthly means on a 1 x 1 degree grid, W m^-2 -- into
%   latitude x longitude x month arrays. Each file is recognised by its own
%   variable_id attribute, so the order of FILES does not matter, and the three
%   must agree on product, version, grid and months or the import refuses.
%
%   WHAT THE GRID MEANS is fixed by the product, not guessed here, and travels
%   in the result (rule R2): the fluxes are defined at a TOA REFERENCE LEVEL
%   20 km above the surface, and the latitudes are GEODETIC (EBAF Ed4.2 Data
%   Quality Summary: the 20 km flux reference level after Loeb et al. 2002, and
%   geodetic zone weights for its global means). cer.referenceHeight_m and
%   cer.latitude carry both.
%
%   'Period' loads only the months whose mid-month time nodes bracket it, one
%   node at or before its start and one at or after its end, which is what a
%   linear interpolation in time over the period needs and nothing more. The
%   whole record is 307 months x 64800 cells x 3 fields: 230 MB as single, so
%   give a period for anything but a climatological study.
%
%   Obtain the files anonymously from ESGF (obs4MIPs, source_id CERES-EBAF-4-2-1),
%   or run setupNonGravForces, which fetches them into <DATAPOOL>/CERES/EBAF and checks each
%   against the SHA-256 ESGF publishes. SYN1deg-3Hour is not read yet (HANDOVER O-3).
%
% INPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   files   1x3   string  -      The rsut, rlut and rsdt obs4MIPs netCDF files, any order
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type      Units  Description
%   ----     ----  ----      -----  -----------
%   Period   1x2   datetime  UTC    [start end], start <= end. Default: every month in the files.
%
% OUTPUT ARGUMENTS
%   Name    Size  Type    Units   Description
%   ----    ----  ----    -----   -----------
%   cer     1x1   struct  -       .product (source_id), .version, .doi, .license,
%                                 .time (1xT datetime, UTC, the files' mid-month nodes),
%                                 .mjd (1xT, the same nodes as MJD UTC),
%                                 .lat_deg (1xNlat centres, ascending), .lon_deg (1xNlon centres,
%                                 ascending, east, a full circle), .latitude "geodetic",
%                                 .referenceHeight_m 20e3,
%                                 .rsut, .rlut, .rsdt (NlatxNlonxT single, W m^-2)
%   aux     1x1   struct  -       .files (1x3, in the order rsut rlut rsdt), .trackingId (1x3),
%                                 .timeIndex (1xT indices into the files' time axis),
%                                 .monthsInFiles, .provenance
%
% ACCURACY
%   Exact: the values are the files' single-precision numbers, only permuted
%   into latitude x longitude x time. The product's own uncertainty governs:
%   monthly regional TOA fluxes carry a few W m^-2 (EBAF DQS), and a monthly mean
%   has no diurnal cycle, which env.earthRadiance restores for the shortwave only,
%   through the instantaneous insolation.
%
% ERRORS
%   NonGravForces:importCERES:inputFiles    not three files, or not one each of rsut, rlut, rsdt
%   NonGravForces:importCERES:dataMissing   a file does not exist, or a loaded value is a fill value,
%                                           NaN or negative
%   NonGravForces:importCERES:dataFormat    not an obs4MIPs CERES EBAF file (a global attribute missing),
%                                           units not W m-2, dimensions not lon/lat/time, time units
%                                           not "days since", non-uniform grid, or the files disagree
%   NonGravForces:importCERES:inputPeriod   Period malformed or not bracketed by the files' months
%
% REFERENCES
%   Loeb N.G. et al. (2018) J. Climate 31(2):895-918 (EBAF Ed4.0); NASA LaRC, CERES_EBAF_Ed4.2
%   Data Quality Summary (2025), sect. on the 20 km reference level and geodetic weights;
%   Loeb N.G., Kato S., Wielicki B.A. (2002) J. Climate 15:3301 (the flux reference level);
%   doi:10.5067/TERRA-AQUA-NOAA20/CERES/EBAF_L3B004.2.1.
%
% EXAMPLE
%   d = fullfile(fileparts(which('setupNonGravForces')), 'data', 'ceres');
%   cer = NonGravForces.io.importCERES(string(fullfile(d, {dir(fullfile(d, '*.nc')).name})));
%   size(cer.rsut), cer.time
%
% See also NonGravForces.env.earthRadiance, NonGravForces.env.radiationSources, setupNonGravForces.
%
% KEYWORDS CERES EBAF albedo outgoing longwave radiation netCDF obs4MIPs Earth radiation pressure
% Created 2026-08-27 by Claude Fable 5 (Anthropic) as a stub; implemented 2026-09-15 by Claude Opus 5 (Anthropic).
arguments
    files string
    options.Period datetime = datetime.empty
end
unit = "importCERES";
need = ["rsut" "rlut" "rsdt"];
files = files(:)';
if numel(files) ~= 3
    NonGravForces.fail(unit, "inputFiles", "got %d file(s) -> pass the three obs4MIPs files rsut, rlut and rsdt of one product", numel(files));
end
missing = files(~isfile(files));
if ~isempty(missing)
    NonGravForces.fail(unit, "dataMissing", "not found: %s -> run setupNonGravForces or fetch the obs4MIPs CERES-EBAF files", strjoin(missing, ", "));
end
vid = strings(1, 3);
for k = 1:3
    vid(k) = localGlobal(unit, files(k), 'variable_id');
end
[tf, order] = ismember(need, vid);
if ~all(tf) || numel(unique(vid)) ~= 3
    NonGravForces.fail(unit, "inputFiles", "the files hold variables %s -> pass exactly one each of rsut, rlut and rsdt", strjoin(vid, ", "));
end
files = files(order);
product = arrayfun(@(f) localGlobal(unit, f, 'source_id'), files);
version = arrayfun(@(f) localGlobal(unit, f, 'source_version_number'), files);
if ~all(startsWith(product, "CERES-EBAF")) || numel(unique(product)) ~= 1 || numel(unique(version)) ~= 1
    NonGravForces.fail(unit, "dataFormat", "source_id %s, version %s -> the three files must be one CERES-EBAF product and version", strjoin(product, "/"), strjoin(version, "/"));
end
% --- axes, identical in the three files ---------------------------------------------------------
lat = ncread(files(1), 'lat')'; lon = ncread(files(1), 'lon')'; tDays = ncread(files(1), 'time')';
tUnits = string(ncreadatt(files(1), 'time', 'units'));
for k = 2:3
    if ~isequal(ncread(files(k), 'lat')', lat) || ~isequal(ncread(files(k), 'lon')', lon) || ~isequal(ncread(files(k), 'time')', tDays) ...
            || string(ncreadatt(files(k), 'time', 'units')) ~= tUnits
        NonGravForces.fail(unit, "dataFormat", "%s and %s differ in their grid or months -> use files of one download", files(1), files(k));
    end
end
% No optional GROUP around the time of day: MATLAB's regexp returned no token from inside a
% (?:...)? group even when it matched, and read "12:30:15" as midnight (measured 2026-09-15).
tok = regexp(tUnits, '^days since (\d{4})-(\d{1,2})-(\d{1,2})[ T]?(\d{0,2}):?(\d{0,2}):?([\d.]*)\s*$', 'tokens', 'once');
if isempty(tok)
    NonGravForces.fail(unit, "dataFormat", "time units '%s' -> expected 'days since yyyy-mm-dd hh:mm:ss' (CF, standard calendar)", tUnits);
end
tok = string(tok);
tok(strlength(tok) == 0) = "0";                                        % a date without a time of day is midnight
t6 = double(tok);
origin = datetime(t6(1), t6(2), t6(3), t6(4), t6(5), t6(6));
dLat = diff(lat); dLon = diff(lon);
if isempty(dLat) || isempty(dLon) || any(abs(dLat - dLat(1)) > 1e-9) || any(abs(dLon - dLon(1)) > 1e-9) || dLat(1) <= 0 || dLon(1) <= 0 ...
        || abs(numel(lon) * dLon(1) - 360) > 1e-9 || lat(1) - dLat(1) / 2 < -90 - 1e-9 || lat(end) + dLat(1) / 2 > 90 + 1e-9
    NonGravForces.fail(unit, "dataFormat", "the grid is not uniform, ascending and global in longitude -> an obs4MIPs 1 x 1 degree EBAF file is expected");
end
nodes = origin + days(tDays);
% --- which months ---------------------------------------------------------------------------------
idx = 1:numel(nodes);
if ~isempty(options.Period)
    p = options.Period;
    if numel(p) ~= 2 || any(isnat(p)) || p(2) < p(1)
        NonGravForces.fail(unit, "inputPeriod", "'Period' must be [start end] with start <= end -> got %d value(s)", numel(p));
    end
    if ~isempty(p.TimeZone), p.TimeZone = 'UTC'; p.TimeZone = ''; end
    k0 = find(nodes <= p(1), 1, 'last'); k1 = find(nodes >= p(2), 1, 'first');
    if isempty(k0) || isempty(k1)
        NonGravForces.fail(unit, "inputPeriod", "the files' months run from %s to %s and cannot bracket %s .. %s -> fetch files that cover the period (interpolation needs a node on each side)", ...
            string(nodes(1)), string(nodes(end)), string(p(1)), string(p(2)));
    end
    idx = k0:k1;
end
% --- values ---------------------------------------------------------------------------------------
cer = struct('product', product(1), 'version', version(1), 'doi', localGlobal(unit, files(1), 'doi'), ...
    'license', localGlobal(unit, files(1), 'license'), 'time', nodes(idx), 'mjd', juliandate(nodes(idx)) - 2400000.5, ...
    'lat_deg', lat, 'lon_deg', lon, 'latitude', "geodetic", 'referenceHeight_m', 20e3, ...
    'rsut', [], 'rlut', [], 'rsdt', []);
tracking = strings(1, 3);
for k = 1:3
    v = need(k);
    u = strtrim(string(ncreadatt(files(k), v, 'units')));
    if u ~= "W m-2"
        NonGravForces.fail(unit, "dataFormat", "%s has units '%s' -> expected 'W m-2'", v, u);
    end
    info = ncinfo(files(k), v);
    dims = string({info.Dimensions.Name});
    iLon = find(dims == "lon"); iLat = find(dims == "lat"); iT = find(dims == "time");
    if numel(dims) ~= 3 || isempty(iLon) || isempty(iLat) || isempty(iT)
        NonGravForces.fail(unit, "dataFormat", "%s has dimensions %s -> expected lon, lat and time", v, strjoin(dims, ", "));
    end
    start = ones(1, 3); count = zeros(1, 3);
    count(iLon) = numel(lon); count(iLat) = numel(lat); count(iT) = numel(idx); start(iT) = idx(1);
    x = ncread(files(k), v, start, count);
    x = permute(x, [iLat iLon iT]);                                    % lat x lon x time
    bad = ~isfinite(x) | x >= 1e19 | x < 0;
    if any(bad(:))
        [~, ~, m] = ind2sub(size(x), find(bad, 1));
        NonGravForces.fail(unit, "dataMissing", "%s holds %d fill, NaN or negative value(s), the first in %s -> choose a period without gaps", ...
            v, nnz(bad), string(nodes(idx(m)), 'yyyy-MM'));
    end
    cer.(v) = single(x);
    tracking(k) = localGlobal(unit, files(k), 'tracking_id');
end
aux = struct('files', files, 'trackingId', tracking, 'timeIndex', idx, 'monthsInFiles', numel(nodes));
if nargout > 1
    aux.provenance = NonGravForces.provenance("io.importCERES", 'Options', struct('files', files, 'period', string(options.Period)));
end
end

function s = localGlobal(unit, file, name)
%LOCALGLOBAL A global attribute as a string, or a dataFormat refusal naming the file.
try
    s = string(ncreadatt(file, '/', name));
catch
    NonGravForces.fail(unit, "dataFormat", "%s has no global attribute %s -> not an obs4MIPs file", file, name);
end
end
