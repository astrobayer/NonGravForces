function [events, aux] = readTHR1B(files, options)
%READTHR1B THR1B thruster-activation reader: timeSeriesPreparation.qc.eventTable impulse priors plus per-thruster on-times.
%
%   [events, aux] = NonGravForces.io.readTHR1B(file)
%   [events, aux] = NonGravForces.io.readTHR1B(files, 'TimeRange', tr)
%
%   Reads GRACE-FO THR1B files through gracefoReading and returns every
%   activation as a row of the canonical event table (timeSeriesPreparation.qc.eventTable):
%   class "impulse", source "THR1B", axis "all", confidence 1, amplitude and
%   sigma NaN (a firing is a certain PRIOR for timeSeriesPreparation.qc.detectImpulses, which
%   estimates the amplitude and assigns the axis from the accelerometer data).
%   tStart is the record epoch in UTC, tEnd = tStart + the longest on-time of
%   the record. Records in which no thruster has an on-time are counted in
%   aux.numRecordsNoFiring and not listed as events.
%
%   WHICH THRUSTER FIRED, AND FOR HOW LONG, is in aux: one column per record,
%   one row per thruster in the order the file header declares (on_time_*
%   variables, milliseconds in the file, seconds here). The thruster-to-axis
%   mapping is deliberately NOT applied: it is in no source the toolbox holds
%   and would be a declared, unmeasured fact (rule R3).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more THR1B files of ONE satellite; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type              Units  Description
%   ----        ----  ----              -----  -----------
%   TimeRange   1x2   double|datetime   -      Forwarded to readL1b. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   events   Nx8   table   -      timeSeriesPreparation.qc.eventTable schema; N records with a firing.
%   aux      1x1   struct  -      .epoch (GPS epoch struct of ALL K records),
%                                 .thruster (1xT string, header names without
%                                 "on_time_"), .onTime_s (TxK), .fired (TxK
%                                 logical), .count (TxK cumulative thrust
%                                 counts, thrust_count_* in header order),
%                                 .accumDur_s (struct .attCtrl 1xK .orbCtrl
%                                 1xK, cumulative), .numRecordsNoFiring,
%                                 .qualflgBits, .product, .satellite, .files,
%                                 .numRecords, .sampling_s, .obsNames, .units,
%                                 .provenance.
%
% ACCURACY
%   Exact. On-times are the file's millisecond integers divided by 1000; the
%   GPS->UTC conversion is referenceSystem.convertTimeScale (18 s in 2020).
%   The event table's datetime resolution (about 1e-7 s near 2020) is coarser
%   than the GPS epoch struct in aux.epoch, which is the one to use for
%   sample-exact work.
%
% ERRORS
%   NonGravForces:readTHR1B:dataProduct     a file is not a THR1B product
%   NonGravForces:readTHR1B:dataSatellite   the files mix satellites
%   NonGravForces:readTHR1B:dataEmpty       no record survives TimeRange
%   NonGravForces:readTHR1B:dataColumn      no on_time_* variable in the header
%   NonGravForces:readTHR1B:dataTimeRef     time_ref is not 'G' (GPS) for every record
%   NonGravForces:readTHR1B:depMissing      timeSeriesPreparation (the event table's owner) not on the path
%   gracefoReading:*, referenceSystem:*     propagated unwrapped
%
% REFERENCES
%   JPL D-56935, THR1B format.
%
% EXAMPLE
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'THR1B_2020-01-01_C_04.txt');
%   [T, aux] = NonGravForces.io.readTHR1B(f); height(T), aux.thruster
%
% See also timeSeriesPreparation.qc.eventTable, timeSeriesPreparation.qc.detectImpulses, NonGravForces.io.readMAS1B.
%
% KEYWORDS THR1B thruster impulse event table prior
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
end
unit = "readTHR1B";
loadOpt = struct('TimeRange', options.TimeRange, 'QualityBitMask', []);
[epoch, obs, raw, aux] = l1bLoad(unit, files, "THR1B", loadOpt);
names = string(raw.obsNames(:)');
onNames = names(startsWith(names, "on_time_"));
if isempty(onNames)
    NonGravForces.fail(unit, "dataColumn", "no on_time_* variable in the header (%s) -> this is not a THR1B product", strjoin(names, " "));
end
if isfield(raw, 'time_ref')
    tr = unique(string(raw.time_ref));
    if ~(isscalar(tr) && tr == "G")
        NonGravForces.fail(unit, "dataTimeRef", "time_ref is %s, expected 'G' (GPS) -> a non-GPS time reference needs a new code path", strjoin(tr, ","));
    end
end
iOn = l1bColumns(unit, raw, onNames);
onTime_s = obs(:, iOn)' / 1000;            % milliseconds in the file
countNames = names(startsWith(names, "thrust_count_"));
iCnt = l1bColumns(unit, raw, countNames, false);
aux.epoch = epoch;
aux.thruster = extractAfter(onNames, "on_time_");
aux.onTime_s = onTime_s;
aux.fired = onTime_s > 0;
aux.count = obs(:, iCnt(iCnt > 0))';
iAcc = l1bColumns(unit, raw, ["accum_dur_att_ctrl","accum_dur_orb_ctrl"], false);
aux.accumDur_s = struct('attCtrl', [], 'orbCtrl', []);
if iAcc(1) > 0, aux.accumDur_s.attCtrl = obs(:, iAcc(1))' / 1000; end
if iAcc(2) > 0, aux.accumDur_s.orbCtrl = obs(:, iAcc(2))' / 1000; end
% Events: records with any on-time, epochs in UTC.
any_fired = any(aux.fired, 1);
aux.numRecordsNoFiring = nnz(~any_fired);
utc = referenceSystem.convertTimeScale(epoch, "UTC");
t = referenceSystem.epochTo(utc, "datetime");
t.TimeZone = 'UTC';
t = t(any_fired(:));
dur = max(onTime_s(:, any_fired), [], 1)';
n = numel(t);
% The event table belongs to timeSeriesPreparation since N-01 (2026-10-14); the deprecated forwarder
% timeSeriesPreparation.qc.eventTable would warn on every THR1B read, so the target is called directly.
if isempty(which('timeSeriesPreparation.qc.eventTable'))
    NonGravForces.fail(unit, "depMissing", "the THR1B events are a timeSeriesPreparation.qc.eventTable and timeSeriesPreparation is not on the path -> run setupTimeSeriesPreparation");
end
events = timeSeriesPreparation.qc.eventTable(t, t + seconds(dur), repmat("impulse", n, 1), repmat("all", n, 1), ...
    nan(n, 1), nan(n, 1), repmat("THR1B", n, 1), ones(n, 1));
aux.provenance = NonGravForces.provenance("io.readTHR1B", 'Files', aux.files, 'Options', options);
end
