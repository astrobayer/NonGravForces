function [epoch, field, aux] = readMAG1B(files, options)
%READMAG1B MAG1B reader: measured magnetic field and magnetorquer currents in the SRF.
%
%   [epoch, field, aux] = NonGravForces.io.readMAG1B(file)
%   [epoch, field, aux] = NonGravForces.io.readMAG1B(files, 'TimeRange', tr)
%
%   Reads GRACE-FO MAG1B files through gracefoReading. Two things travel in this
%   product and they serve opposite purposes: the magnetometer's MEASUREMENT of
%   the ambient field, and the magnetorquers' COMMANDED currents. The second is
%   what makes MAG1B interesting to a non-gravitational force toolbox at all --
%   a current in a coil in the Earth's field is a torque on the satellite, and
%   that torque is not gravitational.
%
%   THE SIX TORQUER COLUMNS ARE THREE AXES, NOT SIX. The file carries an A and a
%   B branch per axis, labelled "positive X" and "negative X" and so on. Measured
%   on the full 2020-01-01 GRACE-C day (172800 records): the A columns are never
%   negative, the B columns are never positive, and A and B are NEVER BOTH
%   NON-ZERO on the same axis -- not once. The branches are already signed, so
%   the net current per axis is their sum, and `field.torquerNet_A` is that sum.
%   The branches are returned untouched beside it, because a reader that only
%   offers a derived quantity has decided something for its caller.
%
%   WHAT THIS READER DOES NOT DO IS TURN CURRENTS INTO A TORQUE. That needs the
%   magnetorquer's turns-area product and its mounting, neither of which is in
%   the file and neither of which the GRACE-FO handbook section this project has
%   read supplies (HANDOVER O-37). The magnetic torque model is WP-23/24 and it
%   will own that conversion; this unit reads.
%
%   THE CALIBRATION COLUMNS ARE RETURNED, NOT APPLIED, and their meaning is not
%   established. On the 2020-01-01 GRACE-C day MF_BCal is exactly [0, 65000,
%   -65000] and torque_cal exactly -137.49 for all 172800 records -- constant,
%   and with values that look like limits or sentinels rather than the scale
%   factors their names suggest. Multiplying a measured field by 65000 because a
%   column is called a calibration factor would be a guess with a clean face.
%   They are in `aux.calibration`, per epoch, exactly as read.
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more MAG1B files of ONE satellite; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name        Size  Type              Units  Description
%   ----        ----  ----              -----  -----------
%   TimeRange   1x2   double|datetime   -      Forwarded to readL1b. Default [].
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   epoch    1x1   struct  -      referenceSystem epoch struct, K entries, "GPS".
%   field    1x1   struct  -      .B_srf_T (3xK, TESLA -- the file is microtesla and this
%                                 converts, because the toolbox is SI throughout)
%                                 .torquerNet_A (3xK, AMPERE, the A + B sum per axis)
%                                 .torquerBranch_A (6xK, ampere, rows 1A 2A 3A 1B 2B 3B,
%                                  as read)
%                                 .magnitude_T (1xK, |B|)
%   aux      1x1   struct  -      .epoch (the epoch struct again), .calibration
%                                 (.magField 3xK, .torquer 1xK, unapplied; NaN when the
%                                 columns are absent), .qualityBits (8xK logical, bit 0
%                                 first), .qualflgBits (Kx8, as the other readers),
%                                 .torquerBranchNames (1x6 string), .product, .satellite,
%                                 .files, .numRecords, .sampling_s, .obsNames, .units,
%                                 .provenance.
%
% ACCURACY
%   Exact (file values) apart from two unit conversions, microtesla to tesla and
%   milliampere to ampere, each a single multiplication by a power of ten.
%
% ERRORS
%   NonGravForces:readMAG1B:dataProduct     a file is not a MAG1B product
%   NonGravForces:readMAG1B:dataSatellite   the files mix satellites
%   NonGravForces:readMAG1B:dataEmpty       no record survives TimeRange
%   NonGravForces:readMAG1B:dataColumn      a MfvX/Y/Z_RAW or torque column is missing
%   NonGravForces:readMAG1B:dataTimeRef     time_ref is not 'G' for every record
%   NonGravForces:readMAG1B:dataRange       |B| above 100 microtesla or a branch current above
%                                           1 A -- an Earth field at 490 km reaches about 50
%                                           microtesla and the measured currents reach 55 mA,
%                                           so either bound means the columns were misread
%   gracefoReading:*                        propagated unwrapped
%
% REFERENCES
%   JPL D-56935 GRACE-FO Level-1 Handbook (MAG1B contents and the SRF);
%   the product's own YAML header, which names every column and its units.
%   HANDOVER.md Part G WP-18 and obligation O-37.
%
% EXAMPLE
%   [ep, f] = NonGravForces.io.readMAG1B("data/samples/MAG1B_2020-01-01_C_04.txt");
%   max(f.magnitude_T) * 1e6            % microtesla, about 50 at GRACE-FO altitude
%   any(f.torquerNet_A ~= 0, 1)         % where the magnetorquers were commanded
%
% See also NonGravForces.io.readACC1B, NonGravForces.io.readTHR1B, NonGravForces.env.magneticField.
%
% KEYWORDS mag1b magnetometer magnetorquer magnetic field torquer current reader
% Created 2026-09-09 by Claude Opus 5 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
end
unit = "readMAG1B";
loadOpt = struct('TimeRange', options.TimeRange, 'QualityBitMask', []);
[epoch, obs, raw, aux] = l1bLoad(unit, files, "MAG1B", loadOpt);
if isfield(raw, 'time_ref')
    tr = unique(string(raw.time_ref));
    if ~(isscalar(tr) && tr == "G")
        NonGravForces.fail(unit, "dataTimeRef", "time_ref is %s, expected 'G' (GPS) -> a non-GPS time reference needs a new code path", strjoin(tr, ","));
    end
end
branchNames = ["torque1A", "torque2A", "torque3A", "torque1B", "torque2B", "torque3B"];
iB = l1bColumns(unit, raw, ["MfvX_RAW", "MfvY_RAW", "MfvZ_RAW"]);
iT = l1bColumns(unit, raw, branchNames);
iC = l1bColumns(unit, raw, ["MF_BCalX", "MF_BCalY", "MF_BCalZ", "torque_cal"], false);

B_uT = obs(:, iB)';                                  % microtesla in the file
branch_mA = obs(:, iT)';                             % milliampere in the file
field = struct();
field.B_srf_T = B_uT * 1e-6;
field.magnitude_T = vecnorm(field.B_srf_T);
field.torquerBranch_A = branch_mA * 1e-3;
% A + B, and the sum is legitimate because the branches never overlap: measured on the full
% 2020-01-01 GRACE-C day, no record has both non-zero on one axis, and the columns are signed.
field.torquerNet_A = field.torquerBranch_A(1:3, :) + field.torquerBranch_A(4:6, :);

if max(field.magnitude_T) > 100e-6
    NonGravForces.fail(unit, "dataRange", "the field magnitude reaches %.1f microtesla -> the Earth's field at 490 km is about 50, so the MfvX/Y/Z columns were misread", 1e6 * max(field.magnitude_T));
end
if max(abs(field.torquerBranch_A), [], 'all') > 1
    NonGravForces.fail(unit, "dataRange", "a magnetorquer branch current reaches %.3f A -> the measured maximum is 0.055 A, so the torque columns were misread", max(abs(field.torquerBranch_A), [], 'all'));
end

aux.epoch = epoch;
aux.torquerBranchNames = branchNames;
K = numel(epoch.sec);
aux.calibration = struct('magField', nan(3, K), 'torquer', nan(1, K));
if all(iC(1:3) > 0), aux.calibration.magField = obs(:, iC(1:3))'; end
if iC(4) > 0, aux.calibration.torquer = obs(:, iC(4))'; end
if isfield(raw, 'qualflgBits')
    aux.qualityBits = logical(raw.qualflgBits)';
else
    aux.qualityBits = false(8, K);
end
% the files read, as every other L1B reader records them (2026-10-31: only their number was recorded)
aux.provenance = NonGravForces.provenance("io.readMAG1B", 'Files', aux.files, 'Options', ...
    struct('files', numel(files), 'records', K));
end
