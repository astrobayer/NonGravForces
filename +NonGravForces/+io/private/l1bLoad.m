function [epoch, obs, raw, common] = l1bLoad(unit, files, expectedProducts, options)
%L1BLOAD Shared loader of the io readers: gracefoReading.readL1bMulti plus the common checks.
%
%   [epoch, obs, raw, common] = l1bLoad(unit, files, expectedProducts, options)
%
%   One code path for every Level-1B reader (rule R1): read, refuse a foreign
%   product, refuse mixed satellites, refuse an empty result, build the GPS
%   epoch struct through referenceSystem, and collect the fields every reader
%   reports. Private to the io package; not a public unit.
%
%   EVERY READ GOES THROUGH readL1bMulti, one file or many (rule R1). It
%   accepts a single file and returns what readL1b returns, so no branch is
%   needed and no aux field is fabricated here. Until gracefoReading 1.0.1
%   the multi-file reader raised MATLAB:badsubscript on an all-empty
%   TimeRange (its issue #6, this repository's O-12), so one file had to go
%   through readL1b to keep dataEmpty reachable; that split is now gone and
%   dataEmpty is reachable for any number of files.
%
%   When TimeRange or the mask removes every record, readL1bMulti first
%   warns gracefoReading:readL1bMulti:noRecordsInRange and then dataEmpty is
%   raised below. The foreign warning is deliberately not suppressed: it
%   names the parsed record count and which filters were active, which the
%   NonGravForces error does not, and foreign diagnostics propagate unwrapped
%   here as foreign errors do.
%
% INPUT ARGUMENTS
%   Name              Size  Type    Units  Description
%   ----              ----  ----    -----  -----------
%   unit              1x1   string  -      Calling unit name, for the error identifier
%   files             1xF   string  -      Files in chronological order
%   expectedProducts  1xP   string  -      Admissible product values, e.g. ["ACC1B","ACT1B","ACH1B"]
%   options           1x1   struct  -      .TimeRange, .QualityBitMask as the caller received them
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, scale "GPS"
%   obs      KxM   double  -      Observation columns as gracefoReading returns them
%   raw      1x1   struct  -      gracefoReading aux with .products and .fileNames (plural) guaranteed
%   common   1x1   struct  -      .product .satellite .files .numRecords .sampling_s
%                                 .obsNames .units .qualflgBits (Kx8 logical)
%
% ERRORS
%   NonGravForces:<unit>:dataProduct     a file is not one of expectedProducts
%   NonGravForces:<unit>:dataSatellite   files of different satellites
%   NonGravForces:<unit>:dataEmpty       no record left after TimeRange / QualityBitMask,
%                                        for any number of files
%   gracefoReading:*, MATLAB:validators:mustBeFile   propagated unwrapped
%
% KEYWORDS private loader L1B
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic); split removed 2026-09-08 by Claude Opus 5 (O-12 discharged).
fwd = {'TimeRange', options.TimeRange, 'QualityBitMask', options.QualityBitMask};
[time, obs, raw] = gracefoReading.readL1bMulti(files, fwd{:});
products = unique(string(raw.products));
bad = products(~ismember(products, expectedProducts));
if ~isempty(bad)
    NonGravForces.fail(unit, "dataProduct", "file(s) carry product %s, expected one of %s -> pass %s files", ...
        strjoin(bad, ","), strjoin(expectedProducts, "|"), strjoin(expectedProducts, "/"));
end
if isempty(time)
    NonGravForces.fail(unit, "dataEmpty", "no record left in %s after TimeRange/QualityBitMask -> widen TimeRange or relax the mask", ...
        strjoin(files, ", "));
end
sat = unique(string(raw.GRACEFO_id));
if numel(sat) ~= 1
    NonGravForces.fail(unit, "dataSatellite", "files mix satellites %s -> read one satellite at a time", strjoin(sat, ","));
end
epoch = referenceSystem.daySecondToEpoch(time, "GPS");
common.product    = products;
common.satellite  = sat;
common.files      = string(files(:)');
common.numRecords = size(time, 1);
if size(time, 1) > 1
    common.sampling_s = median(diff(time(:,1)) * 86400 + diff(time(:,2)));
else
    common.sampling_s = NaN;
end
common.obsNames   = string(raw.obsNames(:)');
common.units      = raw.variableUnits;
if isfield(raw, 'qualflgBits')
    common.qualflgBits = logical(raw.qualflgBits);
else
    common.qualflgBits = false(size(time, 1), 0);
end
end
