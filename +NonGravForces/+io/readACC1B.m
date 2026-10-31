function [epoch, a_srf, aux] = readACC1B(files, options)
%READACC1B ACC1B/ACT1B/ACH1B reader over gracefoReading: 3xK SRF accelerations, GPS epochs, aux.
%
%   [epoch, a_srf, aux] = NonGravForces.io.readACC1B(file)
%   [epoch, a_srf, aux] = NonGravForces.io.readACC1B(files, 'TimeRange', tr, 'QualityBitMask', bits)
%
%   Reads one or more consecutive GRACE-FO accelerometer Level-1B files through
%   gracefoReading.readL1bMulti and returns the linear accelerations in the
%   Science Reference Frame (SRF, which is the toolbox body frame: +X to the
%   K-band horn, +Z to nadir) as a 3xK array, with the epochs as a
%   referenceSystem GPS epoch struct. ACT1B (transplant-calibrated, the
%   GRACE-FO standard), ACH1B (GRACE-D's 1 Hz hybrid transplant from the .ACX
%   archives, ACU + ACS + ACM + AC0, linear accelerations only -- "high rate"
%   until 2026-10-20, which it is not) and classic ACC1B share one layout;
%   the product actually read is reported in aux.product and comes from the
%   file header, never from the file name (rule R3).
%
%   NOTHING IS CALIBRATED HERE. The values are the product's, bias and scale
%   included; calibration is NonGravForces.app.calibrate. Quality flags are
%   returned (aux.qualflgBits), not applied, unless 'QualityBitMask' names the
%   bits whose records are removed (rule R5: the caller decides).
%
% INPUT ARGUMENTS
%   Name     Size  Type    Units  Description
%   ----     ----  ----    -----  -----------
%   files    1xF   string  -      One or more ACC1B/ACT1B/ACH1B files of ONE
%                                 satellite, chronological; .gz accepted.
%
% NAME-VALUE ARGUMENTS
%   Name            Size  Type              Units  Description
%   ----            ----  ----              -----  -----------
%   TimeRange       1x2   double|datetime   -      Inclusive window forwarded
%                                                  to gracefoReading.readL1b:
%                                                  decimal MJD, or datetime on
%                                                  the GPS scale. Default []
%                                                  (whole files).
%   QualityBitMask  1xQ   double            -      qualflg bits (0-based, bit 0
%                                                  rightmost) whose records are
%                                                  removed; forwarded. Default
%                                                  [] (nothing removed).
%
% OUTPUT ARGUMENTS
%   Name     Size  Type    Units   Description
%   ----     ----  ----    -----   -----------
%   epoch    1x1   struct  -       referenceSystem epoch struct, K entries,
%                                  scale "GPS": .day (MJD), .sec (s of day).
%   a_srf    3xK   double  m/s^2   lin_accl_x, _y, _z in the SRF.
%   aux      1x1   struct  -       .angular_rad_s2 (3xK, ang_accl_*),
%                                  .residual_m_s2 (3xK acl_*_res, or [] if the
%                                  product has none), .qualflgBits (Kx8
%                                  logical), .product, .satellite ("C"|"D"),
%                                  .files, .numRecords, .sampling_s (median
%                                  spacing), .obsNames, .units, .provenance.
%
% ACCURACY
%   Exact. Values are the file's doubles (the ASCII carries fewer digits than
%   a double); the epoch split (MJD, seconds of day) is lossless and the time
%   tag round trip is bitwise (tested). No filtering, no interpolation, no
%   gap handling: gaps are reported by gracefoReading.findGaps on the epochs.
%
% ERRORS
%   NonGravForces:readACC1B:dataProduct     a file is not an ACC1B/ACT1B/ACH1B product
%   NonGravForces:readACC1B:dataSatellite   the files mix GRACE-C and GRACE-D
%   NonGravForces:readACC1B:dataEmpty       no record survives TimeRange/QualityBitMask
%   NonGravForces:readACC1B:dataColumn      lin_accl_* or ang_accl_* missing from the header
%   gracefoReading:*                        propagated unwrapped (missing file,
%                                           malformed record, column mismatch)
%
% REFERENCES
%   GRACE-FO Level-1 Data Product User Handbook, JPL D-56935, ACT1B format.
%
% EXAMPLE
%   f = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))), 'data', 'samples', 'ACT1B_2020-01-01_C_04.txt');
%   [ep, a, aux] = NonGravForces.io.readACC1B(f); size(a), aux.product
%
% See also NonGravForces.io.readSCA1B, NonGravForces.io.readGNV1B, gracefoReading.readL1bMulti.
%
% KEYWORDS ACC1B ACT1B accelerometer reader SRF
% Created 2026-09-06 by Claude Fable 5.1 (Anthropic).
arguments
    files (1,:) string {mustBeNonempty}
    options.TimeRange = []
    options.QualityBitMask (1,:) double = []
end
unit = "readACC1B";
[epoch, obs, raw, aux] = l1bLoad(unit, files, ["ACC1B","ACT1B","ACH1B"], options);
iLin = l1bColumns(unit, raw, ["lin_accl_x","lin_accl_y","lin_accl_z"]);
iAng = l1bColumns(unit, raw, ["ang_accl_x","ang_accl_y","ang_accl_z"]);
iRes = l1bColumns(unit, raw, ["acl_x_res","acl_y_res","acl_z_res"], false);
a_srf = obs(:, iLin)';
aux.angular_rad_s2 = obs(:, iAng)';
if all(iRes > 0)
    aux.residual_m_s2 = obs(:, iRes)';
else
    aux.residual_m_s2 = [];
end
aux.provenance = NonGravForces.provenance("io.readACC1B", 'Files', aux.files, 'Options', options);
end
