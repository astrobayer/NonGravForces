function [grid, aux] = importGroopsGrid(file, options)
%IMPORTGROOPSGRID Read a GROOPS binary gridded-data file: global grid, ellipsoid, solid angles and values.
%
%   grid = NonGravForces.io.importGroopsGrid(file)
%   grid = NonGravForces.io.importGroopsGrid(file, 'TimeRange', [mjd0 mjd1])
%   [grid, aux] = NonGravForces.io.importGroopsGrid(___)
%
%   Reads the two container types GROOPS publishes its Earth-radiation data in
%   (https://ftp.tugraz.at/pub/ITSG/groops/data/):
%
%     griddedData            one grid with C value columns. The albedo folder's
%                            earth_ceres_{reflectivity,emissivity}_grid2.5.dat are
%                            144 x 72 with C = 12, a monthly climatology.
%     griddedDataTimeSeries  T+1 time NODES bounding T grids, piecewise constant in
%                            time. earthRadiationPressure/CERES_SYN1deg_*_climatology.dat
%                            are 360 x 180 with 2413 nodes and 2412 monthly grids.
%
%   THE LAYOUT WAS DECODED FROM THE BYTES. GROOPS is GPLv3 and its source is not
%   transcribed here in any language (HANDOVER O-19); no format specification was
%   in hand either, so the structure below is this project's reading of real files
%   and the reference test asserts the three facts that make a misreading visible:
%   the solid angles sum to 4*pi, the value payload fills the file EXACTLY, and for
%   a time series the node count is exactly one more than the grid count.
%
%     'B'  uint64 n  <n bytes: type>  uint64 m  <m bytes: version>  spaces
%     griddedData            : uint64 flag, double a, double b (ellipsoid)
%     griddedDataTimeSeries  : uint64 x 4 (the second is the node count), a, b
%     then four count-prefixed double vectors: longitude and GEOCENTRIC latitude in RADIANS,
%     one radius per latitude row in metres, and one solid angle per point;
%     then uint64 1, double 1.0;
%     time series only: T+1 records of (double fraction, int64 MJD), then uint64;
%     then the values, longitude varying FASTEST within a latitude row.
%
%   SIZE. The SYN1deg files are 1.25 GB each. Without 'TimeRange' every grid is
%   read and the values alone are 1.2 GB in memory; with it, only the grids whose
%   month intersects the range are read, by seeking. Reading two months of SYN1deg
%   costs about 1 MB.
%
% INPUT ARGUMENTS
%   Name    Size  Type    Units  Description
%   ----    ----  ----    -----  -----------
%   file    1x1   string  -      Path to the GROOPS binary file
%
% NAME-VALUE ARGUMENTS
%   Name       Size  Type    Units  Description
%   ----       ----  ----    -----  -----------
%   TimeRange  1x2   double  MJD    Read only the grids whose interval intersects [mjd0 mjd1].
%                                   Time series only; default empty, which reads every grid.
%
% OUTPUT ARGUMENTS
%   Name   Size  Type    Units  Description
%   ----   ----  ----    -----  -----------
%   grid   1x1   struct  -      .type, .version (string); .lon_deg (1xNlon), .lat_deg (1xNlat,
%                               GEOCENTRIC, decreasing -- see below); .radius_m (1xNlat); .area_sr (NlatxNlon);
%                               .values (NlatxNlonxT; for a griddedData T is its C value
%                               columns); .mjd (1x(T+1) node epochs BOUNDING THE
%                               GRIDS RETURNED, so after a 'TimeRange' seek it is the subrange and
%                               not the file's whole node list; [] for griddedData);
%                               .isTimeSeries (logical)
%   aux    1x1   struct  -      .file, .bytes, .ellipsoid (.a, .b), .pointCount, .gridCount,
%                               .gridsRead (1xT index of the grids actually read, into the FILE),
%                               .mjdFile (1x2 first and last node of the file, whatever was read),
%                               .areaSum_sr (must be 4*pi), .provenance
%
% ACCURACY
%   Exact: every value is a float64 read unchanged, and the axes are converted from
%   radians to degrees by one multiplication. Nothing is interpolated here. The
%   solid angles are the file's own and sum to 4*pi to 4e-6 relative (measured
%   2026-09-17 on both 2.5 degree grids), which is the file's discretisation of the
%   sphere, not this reader's error.
%   THE LATITUDES ARE GEOCENTRIC, and this block said "geodetic" until 2026-09-17.
%   The grid is laid out uniformly in GEODETIC latitude and the file stores the
%   GEOCENTRIC latitude of those same points, which is why the stored centres are
%   89.4966, 88.4899, ... rather than 89.5, 88.5, ... and why their spacing runs from
%   0.9933 to 1.0067 instead of being uniform. Two independent measurements settle it:
%     - atan((1 - e^2) tan(89.5 deg)) = 89.49664, the file's first centre, to 1e-7 deg;
%     - .radius_m equals the ellipsoid radius at the stored latitude read as GEOCENTRIC,
%       plus a constant 30000.1 m, to a spread of 0.168 m over all 180 rows. Read as
%       geodetic the same comparison spreads by 71.6 m, i.e. 400 times worse.
%   A consumer that converts its own geocentric latitude to geodetic before looking a
%   value up in this grid misplaces every point by up to 0.19 deg, a fifth of a cell.
%   The longitudes carry single-precision jitter: neighbouring gaps differ by 1.37e-05
%   deg and the first centre is -179.499994 rather than -179.5 (3.2e-08 relative). The
%   axis is uniform in intent, not in its bytes; index it with the average spacing.
%
% ERRORS
%   NonGravForces:importGroopsGrid:inputFile       the file does not exist or cannot be opened
%   NonGravForces:importGroopsGrid:dataFormat      not a GROOPS binary, the solid angles do not
%                                                  match the axes, the payload does not fill the
%                                                  file, or the node, column and grid counts disagree
%   NonGravForces:importGroopsGrid:unsupportedType a GROOPS type this reader does not serve
%   NonGravForces:importGroopsGrid:inputRange      'TimeRange' not two increasing finite MJD,
%                                                  given for a file that carries no time, or
%                                                  covered by no grid of the file
%
% REFERENCES
%   GROOPS data pool, Institute of Geodesy, TU Graz,
%   https://ftp.tugraz.at/pub/ITSG/groops/data/ (albedo and earthRadiationPressure);
%   the CERES products behind them are NASA LaRC's. HANDOVER.md Part G LP-E, WP-E1.
%
% EXAMPLE
%   f = fullfile(getenv('DATAPOOL'), 'EarthRadiationPressure', 'CERES_SYN1deg_albedo_climatology.dat');
%   g = NonGravForces.io.importGroopsGrid(f, 'TimeRange', [58849 58850]);   % 2020-01-01
%   size(g.values)
%
% See also NonGravForces.io.importCERES, NonGravForces.env.earthRadiance, setupNonGravForces.
%
% KEYWORDS groops gridded data ceres albedo longwave earth radiation climatology
% Created 2026-09-17 by Claude Opus 5 (Anthropic).
arguments
    file (1,1) string
    options.TimeRange double = []
end
unit = "importGroopsGrid";
if ~isfile(file)
    NonGravForces.fail(unit, "inputFile", "%s does not exist -> run setupNonGravForces with 'FetchGroops' set to ifMissing and DATAPOOL named", file);
end
tr = options.TimeRange;
if ~isempty(tr) && ~(isequal(size(tr), [1 2]) && all(isfinite(tr)) && tr(2) >= tr(1))
    NonGravForces.fail(unit, "inputRange", "'TimeRange' must be two increasing finite MJD -> got %s", mat2str(tr));
end

fid = fopen(file, 'r', 'l');
if fid < 0
    NonGravForces.fail(unit, "inputFile", "%s cannot be opened -> check the permissions", file);
end
closer = onCleanup(@() fclose(fid));
bytes = dir(file).bytes;

if fread(fid, 1, '*char') ~= 'B'
    NonGravForces.fail(unit, "dataFormat", "%s does not start with the GROOPS binary marker 'B' -> it may be an HTML error page saved as .dat", file);
end
type = localString(fid);
version = localString(fid);
if ~any(type == ["griddedData", "griddedDataTimeSeries"])
    NonGravForces.fail(unit, "unsupportedType", "the file holds GROOPS type '%s' -> this reader serves griddedData and griddedDataTimeSeries", type);
end
isSeries = type == "griddedDataTimeSeries";
localSkipSpaces(fid);

% The preamble differs by type; the node count of a series lives in it.
if isSeries
    pre = fread(fid, 4, 'uint64=>double');
    nNode = pre(2);
else
    fread(fid, 1, 'uint64=>double');
    nNode = 0;
end
ell = fread(fid, 2, 'double');

lon = localVector(fid);
lat = localVector(fid);
rad = localVector(fid);
area = localVector(fid);
nLon = numel(lon); nLat = numel(lat); nPts = numel(area);
if nPts ~= nLon * nLat
    NonGravForces.fail(unit, "dataFormat", "the file has %d longitudes and %d latitudes but %d solid angles -> the geometry does not describe a full grid", nLon, nLat, nPts);
end
fread(fid, 1, 'uint64=>double');                    % a count of 1 ...
fread(fid, 1, 'double');                            % ... and the value 1.0

mjd = [];
if isSeries
    mjd = localNodeEpochs(fid, nNode);
end
% One uint64 stands between the geometry and the values: the COLUMN COUNT in a griddedData
% (12 for the monthly climatologies) and 0 in a time series, where the node count already
% said how many grids follow. It is read either way and cross-checked below.
nCol = fread(fid, 1, 'uint64=>double');
dataStart = ftell(fid);

nGrid = (bytes - dataStart) / 8 / nPts;
if nGrid ~= floor(nGrid) || nGrid < 1
    NonGravForces.fail(unit, "dataFormat", "the values do not fill %s exactly (%g grids of %d points) -> the layout is not what this reader decoded", file, nGrid, nPts);
end
if isSeries && nNode ~= nGrid + 1
    NonGravForces.fail(unit, "dataFormat", "%d time nodes bound %d grids -> a piecewise-constant series needs exactly one more node than grids", nNode, nGrid);
end
if ~isSeries && nCol ~= nGrid
    NonGravForces.fail(unit, "dataFormat", "the file declares %d value columns but %d fit in it -> the layout is not what this reader decoded", nCol, nGrid);
end

want = 1:nGrid;
if ~isempty(tr)
    if ~isSeries
        NonGravForces.fail(unit, "inputRange", "'TimeRange' was given but %s carries no time -> it is a griddedData with %d value columns", file, nGrid);
    end
    want = find(mjd(1:end-1) <= tr(2) & mjd(2:end) >= tr(1));
    if isempty(want)
        NonGravForces.fail(unit, "inputRange", "no grid covers MJD %.4f to %.4f -> the file spans %.4f to %.4f", tr(1), tr(2), mjd(1), mjd(end));
    end
end

values = zeros(nLat, nLon, numel(want));
for k = 1:numel(want)
    fseek(fid, dataStart + (want(k) - 1) * nPts * 8, 'bof');
    v = fread(fid, nPts, 'double');
    values(:, :, k) = reshape(v, nLon, nLat)';      % longitude varies fastest within a latitude row
end

% .mjd IS TRIMMED TO THE GRIDS RETURNED, so the struct says which months it holds. It used to
% carry the file's whole node list beside a seeked subset of values, and a consumer could not tell
% from the struct alone which grid was which -- env.earthRadiance's ceresClimatology backend hit
% exactly that on 2026-09-17 (2413 nodes beside 2 grids). `want` is contiguous by construction, so
% the bounding nodes are one run. The file's full span stays in aux, where it is a fact about the
% FILE rather than about these values (R2: metadata travels with the data it describes).
nodes = mjd; span = [];
if isSeries && ~isempty(mjd)
    nodes = mjd(want(1):want(end) + 1);
    span = mjd([1 end]);
end
grid = struct('type', type, 'version', version, ...
    'lon_deg', rad2deg(lon(:))', 'lat_deg', rad2deg(lat(:))', 'radius_m', rad(:)', ...
    'area_sr', reshape(area, nLon, nLat)', 'values', values, ...
    'mjd', nodes, 'isTimeSeries', isSeries);
aux = struct('file', string(file), 'bytes', bytes, ...
    'ellipsoid', struct('a', ell(1), 'b', ell(2)), 'pointCount', nPts, 'gridCount', nGrid, ...
    'gridsRead', want, 'mjdFile', span, 'areaSum_sr', sum(area), ...
    'provenance', NonGravForces.provenance(unit, 'Files', string(file)));
end

function s = localString(fid)
%LOCALSTRING A length-prefixed string: uint64 count, then that many characters.
n = fread(fid, 1, 'uint64=>double');
s = string(fread(fid, n, '*char')');
end

function localSkipSpaces(fid)
%LOCALSKIPSPACES The header is padded with spaces whose count depends on the type name.
while true
    c = fread(fid, 1, '*char');
    if isempty(c) || c ~= ' '
        fseek(fid, -1, 'cof');
        return
    end
end
end

function v = localVector(fid)
%LOCALVECTOR A count-prefixed vector of doubles.
n = fread(fid, 1, 'uint64=>double');
v = fread(fid, n, 'double');
end

function mjd = localNodeEpochs(fid, nNode)
%LOCALNODEEPOCHS T+1 time nodes, each a double fraction of a day followed by an int64 MJD.
%   Read as two strided streams rather than a loop: 2413 nodes would otherwise be 2413 reads.
%   The file position is left just past the last node, where the caller expects it.
start = ftell(fid);
frac = fread(fid, nNode, 'double', 8);
fseek(fid, start + 8, 'bof');
day = fread(fid, nNode, 'int64=>double', 8);
fseek(fid, start + nNode * 16, 'bof');
mjd = (day + frac)';
end
