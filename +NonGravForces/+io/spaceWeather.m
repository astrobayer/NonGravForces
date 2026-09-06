function sw = spaceWeather(files, options)
%SPACEWEATHER Space-weather provider: CelesTrak SW-All (MSIS drivers) or SET SOLFSMY+DTCFILE (JB2008 drivers), with query(epoch).
%
%   sw = NonGravForces.io.spaceWeather(swAllCsv)                       % family celestrak
%   sw = NonGravForces.io.spaceWeather([solfsmyTxt dtcfileTxt])        % family jb2008
%   d  = sw.query(epoch)                                               % drivers at K epochs
%
%   Loads the solar and geomagnetic indices a thermosphere model needs and
%   returns them with a query function. The two families are never mixed:
%   NRLMSISE-00 / MSIS 2.1 want F10.7 and ap, JB2008 wants F10, S10, M10, Y10
%   and DTC, and the JB2008 indices are not F10.7 substitutes. The family is
%   identified from the FILE CONTENT (header line), never from the file name.
%
%   QUERY. d = sw.query(epoch) takes a referenceSystem epoch struct in any
%   time scale (converted to UTC for the day lookup) and returns 1xK fields:
%
%     celestrak  f107      observed F10.7 of the PREVIOUS UTC day    sfu
%                f107a     81-day centred mean of observed F10.7,
%                          of the day itself                          sfu
%                ap        daily Ap of the day                        2 nT
%                apHistory 7xK: [daily Ap; ap now; ap -3 h; -6 h; -9 h;
%                          mean ap -12..-33 h; mean ap -36..-57 h]    2 nT
%                          (the NRLMSISE-00 AP array, in its order)
%                kp        3-hourly Kp of the epoch, file value / 10  -
%                          (CelesTrak writes thirds as tenths: 0+ = 0.3)
%                dataType  "OBS" | "INT" | "PRD" | "PRM" of the day   -
%     jb2008     f10, f10b        one day earlier                     sfu
%                s10, s10b        one day earlier                     sfu
%                m10, m10b        two days earlier                    sfu
%                y10, y10b        five days earlier                   sfu
%                dtc              at the UTC hour of the epoch         K
%                          (lags as printed in the SOLFSMY header; the
%                          "b" values are the 81-day centred means)
%
%   A query that needs a day outside the table (including the lag and the
%   57 h ap history) raises dataRange; a needed value that is missing (empty
%   field, e.g. Kp in monthly prediction rows) raises dataMissing. The query
%   NEVER returns NaN (HANDOVER E.2 rule 7: a NaN driver silently zeroes or
%   corrupts a density).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One CelesTrak SW-All CSV, or the SET
%                                 SOLFSMY and DTCFILE text files (either order).
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   Format   1x1   string  -      "auto" (default; identified from content),
%                                 "celestrak" or "jb2008" (must then match).
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   sw       1x1   struct  -      .format ("celestrak"|"jb2008"), .mjd (Nx1,
%                                 0h UTC of each day), .span (1x2 MJD),
%                                 .index (struct of Nx1 / Nx8 / Nx24 arrays,
%                                 names below), .files, .query (handle),
%                                 .provenance.
%            celestrak .index:  f107Observed f107Adjusted f107ObservedCentre81
%                                 f107ObservedLast81 apDaily ap3h (Nx8) kp3h
%                                 (Nx8, /10) kpSum sunspotNumber dataType (string)
%            jb2008 .index:     f10 f10b s10 s10b m10 m10b y10 y10b ssrc
%                                 (string) dtc (Nx24, NaN where DTCFILE has no day)
%
% ACCURACY
%   Table values are returned unchanged; nothing is interpolated (daily and
%   3-hourly indices are step functions by definition, DTC is hourly). The
%   81-day means are the published columns, so an excerpt of the file gives
%   the same answers as the whole file.
%
% ERRORS
%   NonGravForces:spaceWeather:inputFiles     wrong number of files for the family
%   NonGravForces:spaceWeather:inputFormat    a file's content matches no known
%                                             family, or not the requested one
%   NonGravForces:spaceWeather:dataColumn     a required column is missing
%   NonGravForces:spaceWeather:dataMalformed  a record does not parse (bad date,
%                                             DTC row without 24 values, JD and
%                                             YYYY/DDD disagree)
%   NonGravForces:spaceWeather:inputEpoch     (query) not an epoch struct
%   NonGravForces:spaceWeather:dataRange      (query) a needed day is not in the table
%   NonGravForces:spaceWeather:dataMissing    (query) a needed value is empty
%   MATLAB:validators:mustBeFile              a file does not exist
%   referenceSystem:*                         propagated unwrapped
%
% REFERENCES
%   CelesTrak SW-All format, https://celestrak.org/SpaceData/ ; Picone et al.
%   (2002) NRLMSISE-00 (AP array, F107 previous day, F107A 81-day centred);
%   Bowman et al. (2008) AIAA 2008-6438, JB2008 driver header (index lags);
%   SET SOLFSMY.TXT header (lags restated).
%
% EXAMPLE
%   root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
%   sw = NonGravForces.io.spaceWeather(fullfile(root, 'data', 'spaceweather', 'SW-All_2019-09_2020-03.csv'));
%   d = sw.query(referenceSystem.mjdToEpoch(58849.5, "UTC"));   % d.f107 = 71.7 (2019-12-31)
%
% See also NonGravForces.env.density, referenceSystem.mjdToEpoch.
%
% KEYWORDS space weather F10.7 Ap Kp S10 M10 Y10 DTC JB2008 MSIS
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty, mustBeFile}
    options.Format (1,1) string {mustBeMember(options.Format, ["auto","celestrak","jb2008"])} = "auto"
end
unit = "spaceWeather";
kinds = arrayfun(@(f) localKind(f, unit), files);
if options.Format == "auto"
    if isscalar(files) && kinds == "celestrak"
        fam = "celestrak";
    elseif numel(files) == 2 && all(ismember(["solfsmy","dtc"], kinds))
        fam = "jb2008";
    else
        NonGravForces.fail(unit, "inputFiles", "got %d file(s) of kind %s -> pass one SW-All.csv (celestrak) or SOLFSMY.TXT + DTCFILE.TXT (jb2008)", ...
            numel(files), strjoin(kinds, ","));
    end
else
    fam = options.Format;
    want = struct('celestrak', "celestrak", 'jb2008', ["solfsmy","dtc"]);
    if ~(numel(files) == numel(want.(fam)) && all(ismember(want.(fam), kinds)))
        NonGravForces.fail(unit, "inputFormat", "Format %s needs file kind(s) %s, got %s -> pass the matching files or use Format auto", ...
            fam, strjoin(want.(fam), "+"), strjoin(kinds, ","));
    end
end
sw.format = fam;
sw.files = files;
switch fam
    case "celestrak"
        [sw.mjd, sw.index] = localReadCelestrak(files, unit);
        data = struct('mjd', sw.mjd, 'index', sw.index);
        sw.query = @(epoch) localQueryCelestrak(epoch, data, unit);
    case "jb2008"
        [sw.mjd, sw.index] = localReadJb2008(files(kinds == "solfsmy"), files(kinds == "dtc"), unit);
        data = struct('mjd', sw.mjd, 'index', sw.index);
        sw.query = @(epoch) localQueryJb2008(epoch, data, unit);
end
sw.span = [sw.mjd(1) sw.mjd(end)];
sw.provenance = NonGravForces.provenance("io.spaceWeather", 'Files', files, 'Options', options);
end

% ------------------------------------------------------------------------
function kind = localKind(file, unit)
% The family is read from the first lines of the file, never from its name.
fid = fopen(file, 'r');
if fid < 0
    NonGravForces.fail(unit, "inputFormat", "%s cannot be opened -> check the path and permissions", file);
end
head = strings(0, 1);
while numel(head) < 6
    l = fgetl(fid);
    if ~ischar(l), break; end
    if ~isempty(strtrim(l)), head(end+1, 1) = string(l); end %#ok<AGROW>
end
fclose(fid);
if ~isempty(head) && startsWith(head(1), "DATE,BSRN,ND,")
    kind = "celestrak";
elseif any(startsWith(head, "DTC "))
    kind = "dtc";
elseif any(contains(head, "YYYY DDD") & startsWith(head, "#"))
    kind = "solfsmy";
else
    NonGravForces.fail(unit, "inputFormat", "%s is neither CelesTrak SW-All (header DATE,BSRN,ND,...), SET SOLFSMY (# YYYY DDD ...) nor SET DTCFILE (DTC yyyy ddd ...) -> pass one of those", file);
end
end

function [mjd, index] = localReadCelestrak(file, unit)
T = readtable(file, 'Delimiter', ',', 'TextType', 'string', 'VariableNamingRule', 'preserve');
need = ["DATE","KP1","KP2","KP3","KP4","KP5","KP6","KP7","KP8","KP_SUM","AP1","AP2","AP3","AP4","AP5","AP6","AP7","AP8","AP_AVG","ISN","F10.7_OBS","F10.7_ADJ","F10.7_DATA_TYPE","F10.7_OBS_CENTER81","F10.7_OBS_LAST81"];
have = string(T.Properties.VariableNames);
miss = need(~ismember(need, have));
if ~isempty(miss)
    NonGravForces.fail(unit, "dataColumn", "SW-All columns %s missing (have %s) -> use the CelesTrak SW-All.csv layout", strjoin(miss, ","), strjoin(have, ","));
end
d = T.DATE;
if ~isdatetime(d)
    d = datetime(string(d), 'InputFormat', 'yyyy-MM-dd');
end
if any(isnat(d))
    NonGravForces.fail(unit, "dataMalformed", "row %d of %s has no parsable DATE -> repair or re-download the file", find(isnat(d), 1), file);
end
d.TimeZone = '';
ep = referenceSystem.datetimeToEpoch(d, "UTC");
mjd = double(ep.day(:));
if any(diff(mjd) <= 0)
    NonGravForces.fail(unit, "dataMalformed", "dates in %s are not strictly increasing at row %d -> the file is corrupt", file, find(diff(mjd) <= 0, 1) + 1);
end
num = @(n) localNumeric(T.(n), n, unit);
index.f107Observed         = num("F10.7_OBS");
index.f107Adjusted         = num("F10.7_ADJ");
index.f107ObservedCentre81 = num("F10.7_OBS_CENTER81");
index.f107ObservedLast81   = num("F10.7_OBS_LAST81");
index.apDaily              = num("AP_AVG");
index.ap3h = nan(numel(mjd), 8); index.kp3h = nan(numel(mjd), 8);
for k = 1:8
    index.ap3h(:, k) = num("AP" + k);
    index.kp3h(:, k) = num("KP" + k) / 10;
end
index.kpSum         = num("KP_SUM") / 10;
index.sunspotNumber = num("ISN");
index.dataType      = string(T.("F10.7_DATA_TYPE"));
end

function x = localNumeric(col, name, unit)
if isnumeric(col)
    x = double(col(:));
else
    x = double(str2double(string(col(:))));   % empty fields -> NaN
end
if ~isvector(x) && ~isempty(x)
    NonGravForces.fail(unit, "dataColumn", "column %s is not a vector -> check the file layout", name);
end
end

function [mjd, index] = localReadJb2008(solfsmy, dtcfile, unit)
% SOLFSMY: comment lines '#', then YYYY DDD JulianDay F10 F81c S10 S81c M10 M81c Y10 Y81c Ssrc
lines = localDataLines(solfsmy, "#");
tok = regexp(lines, '^\s*(\d{4})\s+(\d{1,3})\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s*$', 'tokens', 'once');
bad = find(cellfun(@isempty, tok), 1);
if ~isempty(bad)
    NonGravForces.fail(unit, "dataMalformed", "SOLFSMY data line %d does not have YYYY DDD JD and 8 values plus Ssrc -> check the file: %s", bad, lines(bad));
end
M = vertcat(tok{:});
yr = str2double(M(:, 1)); doy = str2double(M(:, 2)); jd = str2double(M(:, 3));
vals = str2double(M(:, 4:11));
mjd = localMjdFromDoy(yr, doy);
if any(abs(jd - 2400000.5 - 0.5 - mjd) > 1e-6)
    k = find(abs(jd - 2400000.5 - 0.5 - mjd) > 1e-6, 1);
    NonGravForces.fail(unit, "dataMalformed", "SOLFSMY line %d: JulianDay %.1f disagrees with YYYY DDD %d %d -> the file is corrupt", k, jd(k), yr(k), doy(k));
end
if any(diff(mjd) <= 0)
    NonGravForces.fail(unit, "dataMalformed", "SOLFSMY days are not strictly increasing at line %d -> the file is corrupt", find(diff(mjd) <= 0, 1) + 1);
end
names = ["f10","f10b","s10","s10b","m10","m10b","y10","y10b"];
for k = 1:8, index.(names(k)) = vals(:, k); end
index.ssrc = M(:, 12);
% DTCFILE: DTC YYYY DDD v1 .. v24 (hourly)
dl = localDataLines(dtcfile, "#");
tokD = regexp(dl, '^\s*DTC\s+(\d{4})\s+(\d{1,3})((?:\s+\S+){24})\s*$', 'tokens', 'once');
bad = find(cellfun(@isempty, tokD), 1);
if ~isempty(bad)
    NonGravForces.fail(unit, "dataMalformed", "DTCFILE line %d is not 'DTC YYYY DDD' followed by 24 values -> check the file: %s", bad, dl(bad));
end
D = vertcat(tokD{:});
mjdD = localMjdFromDoy(str2double(D(:, 1)), str2double(D(:, 2)));
dtcVals = cell2mat(cellfun(@(s) str2double(split(strtrim(s)))', D(:, 3), 'UniformOutput', false));
index.dtc = nan(numel(mjd), 24);
[tf, loc] = ismember(mjd, mjdD);
index.dtc(tf, :) = dtcVals(loc(tf), :);
end

function lines = localDataLines(file, comment)
txt = string(fileread(file));
lines = splitlines(txt);
lines = strtrim(lines);
lines = lines(lines ~= "" & ~startsWith(lines, comment));
end

function mjd = localMjdFromDoy(yr, doy)
ep = referenceSystem.datetimeToEpoch(datetime(yr, 1, 1) + caldays(doy - 1), "UTC");
mjd = double(ep.day(:));
end

% ------------------------------------------------------------------------
function [dayQ, secQ] = localQueryDay(epoch, unit)
if ~(isstruct(epoch) && all(isfield(epoch, ["day","sec","scale"])))
    NonGravForces.fail(unit, "inputEpoch", "query takes a referenceSystem epoch struct (.day .sec .scale) -> build one with mjdToEpoch or daySecondToEpoch");
end
if string(epoch.scale) ~= "UTC"
    epoch = referenceSystem.convertTimeScale(epoch, "UTC");
end
dayQ = double(epoch.day(:)'); secQ = double(epoch.sec(:)');
end

function rows = localRows(mjd, days, unit, what)
[tf, rows] = ismember(days, mjd);
if ~all(tf(:))
    k = find(~tf, 1);
    NonGravForces.fail(unit, "dataRange", "%s needs MJD %d (%s), the table covers %d..%d -> load a file spanning the arc plus its lags", ...
        what, days(k), string(datetime(days(k), 'ConvertFrom', 'modifiedjuliandate', 'Format', 'yyyy-MM-dd')), mjd(1), mjd(end));
end
end

function localNoNaN(x, name, days, unit)
if any(isnan(x(:)))
    [~, k] = find(isnan(x), 1);
    NonGravForces.fail(unit, "dataMissing", "%s is missing in the table for MJD %d (%s) -> use a file with observed values for that day", ...
        name, days(k), string(datetime(days(k), 'ConvertFrom', 'modifiedjuliandate', 'Format', 'yyyy-MM-dd')));
end
end

function d = localQueryCelestrak(epoch, data, unit)
[dayQ, secQ] = localQueryDay(epoch, unit);
idx = data.index;
r0 = localRows(data.mjd, dayQ, unit, "celestrak query");
r1 = localRows(data.mjd, dayQ - 1, unit, "previous-day F10.7");
slot = floor(secQ / 10800) + 1;                    % 1..8, the 3-h interval of the epoch
K = numel(dayQ);
% 32 consecutive 3-h ap values ending with the day of the epoch: days d-3..d
rHist = zeros(4, K);
for j = 1:4, rHist(j, :) = localRows(data.mjd, dayQ - 4 + j, unit, "ap history (57 h)"); end
A = zeros(32, K);
for j = 1:4, A(8*(j-1)+(1:8), :) = idx.ap3h(rHist(j, :), :)'; end
cur = 24 + slot;                                    % index of the current interval in A
pick = @(off) A(sub2ind(size(A), cur - off, 1:K));
apNow = pick(0); ap3 = pick(1); ap6 = pick(2); ap9 = pick(3);
m1 = zeros(1, K); m2 = zeros(1, K);
for o = 4:11,  m1 = m1 + pick(o); end
for o = 12:19, m2 = m2 + pick(o); end
d.f107  = idx.f107Observed(r1)';
d.f107a = idx.f107ObservedCentre81(r0)';
d.ap    = idx.apDaily(r0)';
d.apHistory = [d.ap; apNow; ap3; ap6; ap9; m1 / 8; m2 / 8];
d.kp    = idx.kp3h(sub2ind(size(idx.kp3h), r0, slot));
d.dataType = idx.dataType(r0)';
localNoNaN(d.f107, "F10.7_OBS (previous day)", dayQ - 1, unit);
localNoNaN(d.f107a, "F10.7_OBS_CENTER81", dayQ, unit);
localNoNaN(d.ap, "AP_AVG", dayQ, unit);
localNoNaN(d.apHistory, "3-hourly ap (57 h history)", dayQ, unit);
localNoNaN(d.kp, "3-hourly Kp", dayQ, unit);
end

function d = localQueryJb2008(epoch, data, unit)
[dayQ, secQ] = localQueryDay(epoch, unit);
idx = data.index;
r0 = localRows(data.mjd, dayQ, unit, "jb2008 query");
r1 = localRows(data.mjd, dayQ - 1, unit, "F10/S10 lag (1 d)");
r2 = localRows(data.mjd, dayQ - 2, unit, "M10 lag (2 d)");
r5 = localRows(data.mjd, dayQ - 5, unit, "Y10 lag (5 d)");
hour = floor(secQ / 3600) + 1;                      % 1..24
d.f10  = idx.f10(r1)';  d.f10b = idx.f10b(r1)';
d.s10  = idx.s10(r1)';  d.s10b = idx.s10b(r1)';
d.m10  = idx.m10(r2)';  d.m10b = idx.m10b(r2)';
d.y10  = idx.y10(r5)';  d.y10b = idx.y10b(r5)';
d.dtc  = idx.dtc(sub2ind(size(idx.dtc), r0, hour));
localNoNaN([d.f10; d.f10b; d.s10; d.s10b], "F10/S10 (1 d lag)", dayQ - 1, unit);
localNoNaN([d.m10; d.m10b], "M10 (2 d lag)", dayQ - 2, unit);
localNoNaN([d.y10; d.y10b], "Y10 (5 d lag)", dayQ - 5, unit);
localNoNaN(d.dtc, "DTC (hourly)", dayQ, unit);
end
