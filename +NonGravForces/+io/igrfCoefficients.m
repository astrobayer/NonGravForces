function [coeff, aux] = igrfCoefficients(file, options)
%IGRFCOEFFICIENTS IGRF-14 Schmidt semi-normalised spherical-harmonic coefficients, with the model's own time rule.
%
%   coeff = NonGravForces.io.igrfCoefficients(file)
%   coeff = NonGravForces.io.igrfCoefficients(file, 'Epoch', ep)
%   coeff = NonGravForces.io.igrfCoefficients(file, 'Epoch', 2020.5)
%
%   Reads the IAGA V-MOD coefficient file (`igrf14coeffs.txt`) into arrays a
%   spherical-harmonic evaluation can use directly. The file publishes the main
%   field at 26 epochs from 1900.0 to 2025.0 in five-year steps, degree 1 to 13,
%   plus a secular variation for 2025-2030 to degree 8, all in nanotesla.
%
%   THE TIME RULE IS THE MODEL'S, NOT THIS READER'S INVENTION, and it is applied
%   here rather than left to every caller because getting it wrong is silent:
%   linear interpolation between the two bracketing epochs, and beyond the last
%   epoch the secular variation extrapolates, g(t) = g(t_last) + SV (t - t_last).
%   That extrapolation is valid for five years and this refuses past it. Before
%   1900.0 it refuses too: the model does not reach there and returning the
%   1900 field for 1850 would be an answer with no basis.
%
%   COEFFICIENTS ARE INDEXED [n+1, m+1], so g(2,1) is g_1^0 and g(2,2) is g_1^1.
%   The zeroth row and the upper triangle are zero and stay in the array: a
%   degree-ordered loop reads them as zero and needs no special case, which is
%   worth more than the memory it costs.
%
%   UNITS STAY NANOTESLA HERE. The rest of the toolbox is SI and
%   NonGravForces.env.magneticField returns tesla, but the coefficients are
%   published in nanotesla and a reader that silently rescales its source makes
%   every published value impossible to check against the file. The conversion
%   belongs to the evaluation, and it happens once.
%
% INPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   file    1x1   string  -      Path to igrf14coeffs.txt (IAGA V-MOD format).
%
% NAME-VALUE ARGUMENTS
%   Name     Size  Type            Units  Description
%   ----     ----  ----            -----  -----------
%   Epoch    1x1   struct|double   -      A referenceSystem epoch struct or a decimal year.
%                                         Given, coeff.g and coeff.h are the interpolated
%                                         field at that time; absent, only the tables are
%                                         returned and coeff.g is empty.
%
% OUTPUT ARGUMENTS
%   Name    Size  Type    Units    Description
%   ----    ----  ----    -----    -----------
%   coeff   1x1   struct  nT       .gTable, .hTable ((nMax+1)x(nMax+1)xE, one page per epoch)
%                                  .epochs (1xE decimal years, 1900.0 ... 2025.0)
%                                  .svG, .svH ((nMax+1)x(nMax+1), nT/year, degree <= svMax)
%                                  .nMax (13), .svMax (8)
%                                  .g, .h ((nMax+1)x(nMax+1) at 'Epoch', else [])
%                                  .year (decimal year used, else NaN)
%   aux     1x1   struct  -        .file, .generation, .numCoefficients, .validFrom, .validTo,
%                                  .header (the two comment lines), .provenance
%
% ACCURACY
%   Exact: the parsed values are the file's, and the interpolation is one linear
%   blend of two published numbers. Nothing is rounded, rescaled or reordered.
%
% ERRORS
%   NonGravForces:igrfCoefficients:dataMissing   the file does not exist
%   NonGravForces:igrfCoefficients:dataFormat    no "g/h n m ..." header line, or a row whose
%                                                field count disagrees with that header
%   NonGravForces:igrfCoefficients:dataColumn    a row is neither g nor h, or n and m are not
%                                                a valid Schmidt pair (m <= n, n >= 1)
%   NonGravForces:igrfCoefficients:inputRange    'Epoch' before the first epoch or more than
%                                                five years past the last
%   NonGravForces:igrfCoefficients:numericNaN    a coefficient does not parse as a number
%
% REFERENCES
%   Alken P. et al. (2021) Earth Planets Space 73:49 (IGRF-13 description; the
%   IGRF-14 release follows the same conventions); IAGA Division V-MOD,
%   https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt (the file this
%   reads, byte-identical to the copy recorded in data/igrf/README.md);
%   Winch D. et al. (2005) Geophys. J. Int. 160:487 (Schmidt semi-normalisation).
%   HANDOVER.md Part G WP-23.
%
% EXAMPLE
%   c = NonGravForces.io.igrfCoefficients("data/igrf/igrf14coeffs.txt", 'Epoch', 2020.0);
%   c.g(2, 1)                      % g_1^0 = -29403.41 nT, the axial dipole
%   c.svG(2, 1)                    % 12.6 nT/year: the dipole is weakening
%
% See also NonGravForces.env.magneticField, NonGravForces.io.readMAG1B.
%
% KEYWORDS igrf geomagnetic coefficients spherical harmonics schmidt secular variation
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    file (1,1) string
    options.Epoch = []
end
unit = "igrfCoefficients";
if ~isfile(file)
    NonGravForces.fail(unit, "dataMissing", "the IGRF coefficient file %s does not exist -> setupNonGravForces fetches it, or download it from https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf14coeffs.txt", file);
end
lines = splitlines(string(fileread(file)));
lines = lines(strlength(strtrim(lines)) > 0);
isComment = startsWith(strtrim(lines), "#");
header = lines(isComment);
body = lines(~isComment);
% The "g/h n m 1900.0 ..." line carries the epochs; the "c/s deg ord ..." line above it is a
% label row and is skipped. Finding the epochs in the file rather than assuming 1900:5:2025
% means a future generation with another span parses without a code change.
iHdr = find(startsWith(strtrim(body), "g/h"), 1);
if isempty(iHdr)
    NonGravForces.fail(unit, "dataFormat", "no 'g/h n m ...' header line in %s -> this is not an IAGA V-MOD coefficient file", file);
end
hdr = split(strtrim(body(iHdr)));
nCol = numel(hdr);
epochs = str2double(hdr(4:end-1))';            % the last column is the secular variation
if any(isnan(epochs)) || isempty(epochs)
    NonGravForces.fail(unit, "dataFormat", "the epoch columns of the header do not parse as years -> got %s", strjoin(hdr(4:end-1), " "));
end
rows = body(iHdr+1:end);
K = numel(rows);
gh = strings(K, 1); nDeg = zeros(K, 1); mOrd = zeros(K, 1); vals = zeros(K, numel(epochs) + 1);
for k = 1:K
    p = split(strtrim(rows(k)));
    if numel(p) ~= nCol
        NonGravForces.fail(unit, "dataFormat", "row %d has %d fields, the header declares %d -> the file is truncated or a column is missing", k, numel(p), nCol);
    end
    gh(k) = p(1);
    nDeg(k) = str2double(p(2));
    mOrd(k) = str2double(p(3));
    v = str2double(p(4:end));
    if any(isnan(v))
        NonGravForces.fail(unit, "numericNaN", "row %d holds a value that is not a number -> %s", k, strjoin(p(4:end), " "));
    end
    vals(k, :) = v';
end
if ~all(gh == "g" | gh == "h")
    NonGravForces.fail(unit, "dataColumn", "a row is labelled neither g nor h -> the first column must say which harmonic it is");
end
if any(nDeg < 1) || any(mOrd < 0) || any(mOrd > nDeg)
    NonGravForces.fail(unit, "dataColumn", "a row has an invalid degree/order pair -> Schmidt coefficients need n >= 1 and 0 <= m <= n");
end

nMax = max(nDeg);
E = numel(epochs);
coeff = struct();
coeff.gTable = zeros(nMax + 1, nMax + 1, E);
coeff.hTable = zeros(nMax + 1, nMax + 1, E);
coeff.svG = zeros(nMax + 1, nMax + 1);
coeff.svH = zeros(nMax + 1, nMax + 1);
for k = 1:K
    i = nDeg(k) + 1; j = mOrd(k) + 1;
    if gh(k) == "g"
        coeff.gTable(i, j, :) = reshape(vals(k, 1:E), 1, 1, E);
        coeff.svG(i, j) = vals(k, end);
    else
        coeff.hTable(i, j, :) = reshape(vals(k, 1:E), 1, 1, E);
        coeff.svH(i, j) = vals(k, end);
    end
end
coeff.epochs = epochs;
coeff.nMax = nMax;
svRows = nDeg(vals(:, end) ~= 0);
coeff.svMax = max([svRows; 0]);
coeff.g = [];
coeff.h = [];
coeff.year = NaN;

if ~isempty(options.Epoch)
    y = localDecimalYear(unit, options.Epoch);
    if y < epochs(1)
        NonGravForces.fail(unit, "inputRange", "epoch %.4f is before the model's first epoch %.1f -> IGRF does not reach there, and returning the earliest field instead would be an answer with no basis", y, epochs(1));
    end
    if y > epochs(end) + 5
        NonGravForces.fail(unit, "inputRange", "epoch %.4f is more than five years past the last epoch %.1f -> the secular variation is published for %.0f-%.0f only", y, epochs(end), epochs(end), epochs(end) + 5);
    end
    if y >= epochs(end)
        dt = y - epochs(end);
        coeff.g = coeff.gTable(:, :, end) + coeff.svG * dt;
        coeff.h = coeff.hTable(:, :, end) + coeff.svH * dt;
    else
        i2 = find(epochs > y, 1);
        i1 = i2 - 1;
        w = (y - epochs(i1)) / (epochs(i2) - epochs(i1));
        coeff.g = (1 - w) * coeff.gTable(:, :, i1) + w * coeff.gTable(:, :, i2);
        coeff.h = (1 - w) * coeff.hTable(:, :, i1) + w * coeff.hTable(:, :, i2);
    end
    coeff.year = y;
end

gen = regexp(strjoin(header, " "), '(\d+)(?:th|st|nd|rd) Generation', 'tokens', 'once');
aux = struct('file', string(file), 'generation', NaN, 'numCoefficients', K, ...
    'validFrom', epochs(1), 'validTo', epochs(end) + 5, 'header', header(:)', ...
    'provenance', NonGravForces.provenance("io.igrfCoefficients", 'Options', ...
        struct('nMax', nMax, 'epochs', E, 'year', coeff.year)));
if ~isempty(gen), aux.generation = str2double(gen{1}); end
end

function y = localDecimalYear(unit, ep)
%LOCALDECIMALYEAR A referenceSystem epoch or a plain decimal year, both accepted.
% The conversion goes through referenceSystem rather than a days-per-year constant, because
% "decimal year" is only unambiguous once someone has decided what a leap year does to it.
if isnumeric(ep)
    if ~isscalar(ep) || ~isfinite(ep)
        NonGravForces.fail(unit, "inputRange", "'Epoch' as a number must be one finite decimal year -> got %s", mat2str(size(ep)));
    end
    y = ep;
    return
end
t = referenceSystem.epochTo(referenceSystem.convertTimeScale(ep, "UTC"), "datetime");
t.TimeZone = 'UTC';
y0 = datetime(year(t), 1, 1, 'TimeZone', 'UTC');
y1 = datetime(year(t) + 1, 1, 1, 'TimeZone', 'UTC');
y = year(t) + days(t - y0) / days(y1 - y0);
end
